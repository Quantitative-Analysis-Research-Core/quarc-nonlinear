function [tau, curve, info] = ami_kde(x, maxlag, opts)
%AMI_KDE Average mutual information versus lag, kernel density estimator.
%   TAU = AMI_KDE(X,maxlag) returns the lag of the first local minimum of the
%   average mutual information of X with itself over lags 0:maxlag, using Gaussian
%   kernel density estimates of the marginal and joint densities.
%
%   [TAU,CURVE] = AMI_KDE(X,maxlag) also returns the AMI curve as an (maxlag+1)-by-2
%   array of [lag, ami], in bits.
%
%   [TAU,CURVE,INFO] = AMI_KDE(X,maxlag) also returns a struct with fields
%   allMinima, usedFallback, fractionLag, bandwidth and estimator.
%
%   ___ = AMI_KDE(X,maxlag,fraction=F) sets the fallback threshold, as a fraction
%   of AMI at lag 0, used when the curve has no local minimum. Default 0.2.
%
%   ___ = AMI_KDE(X,maxlag,bandwidth="adaptive") lets each kernel's width
%   vary with the local density (Abramson 1982): narrow where points are
%   dense, wide in the tails. Default "fixed". About twice the cost.
%
%   ___ = AMI_KDE(X,maxlag,chunk=C) sets how many rows of the pairwise kernel
%   are evaluated at once. Default 512. Lower it to reduce peak memory; it
%   does not change the result.
%
%   Input Arguments
%      X  time series, real column vector, no NaN
%      maxlag  maximum lag, positive integer, less than numel(X)
%
%   Notes
%   Marginals use a univariate Gaussian kernel with bandwidth std(X)/N^(1/6).
%   The joint density uses a bivariate Gaussian kernel whose covariance
%   carries the sample correlation. Densities are evaluated at the data
%   points and AMI is the mean log2 ratio.
%
%   With bandwidth="adaptive", a fixed-bandwidth estimate is used as a pilot
%   and kernel i is scaled by lambda_i = sqrt(g / pilot_i), where g is the
%   geometric mean of the pilot. The joint kernel keeps its correlation and
%   scales both axes by the same lambda_i. A fixed bandwidth oversmooths
%   the dense core of a skewed or multimodal distribution and undersmooths
%   its tails; the adaptive one corrects both, which matters for chaotic
%   series whose densities are far from Gaussian.
%
%   Cost is O(maxlag*N^2) kernel evaluations, substantially slower than
%   AMI_HISTOGRAM, in exchange for lower bias. Peak memory is O(ChunkSize*N).
%
%   Runs on base MATLAB. No toolbox required.
%
%   Examples
%      tau = ami_kde(x, 50);
%      tau = ami_kde(x, 50, chunk=128);   % lower peak memory
%      tau = ami_kde(x, 50, bandwidth="adaptive");
%
%   References
%      Abramson, I. S. (1982). On bandwidth variation in kernel estimates --
%      a square root law. The Annals of Statistics, 10(4), 1217-1223.
%
%      Thomas, R. D., Moses, N. C., Semple, E. A. and Strang, A. J. (2014).
%      An efficient algorithm for the computation of average mutual
%      information. Behavior Research Methods.
%
%      Fraser, A. M. and Swinney, H. maxlag. (1986). Independent coordinates for
%      strange attractors from mutual information. Physical Review A, 33(2),
%      1134-1140.
%
%   See also AMI, AMI_HISTOGRAM.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    x  (:,1) double {mustBeNonempty}
    maxlag  (1,1) double {mustBePositive, mustBeInteger}
    opts.fraction  (1,1) double {mustBePositive} = 0.2
    opts.chunk (1,1) double {mustBePositive, mustBeInteger} = 512
    opts.bandwidth (1,1) string {mustBeMember(opts.bandwidth, ["fixed", "adaptive"])} = "fixed"
end

adaptive = opts.bandwidth == "adaptive";

if anynan(x)
    error('ami_kde:nanInput', ...
          'x contains NaN. Mutual information is undefined for missing data.');
end
if maxlag >= numel(x)
    error('ami_kde:lagTooLarge', ...
          'L (%d) must be smaller than the number of samples (%d).', maxlag, numel(x));
end

curve = zeros(maxlag + 1, 2);
curve(1, :) = [0, mutualInformationBits(x, x, opts.chunk, adaptive)];
for lag = 1:maxlag
    a = x(1:end-lag);
    b = x(lag+1:end);
    curve(lag + 1, :) = [lag, mutualInformationBits(a, b, opts.chunk, adaptive)];
end

[tau, info] = ami_first_minimum(curve, opts.fraction, "ami_kde");
info.bandwidth = opts.bandwidth;
if adaptive
    info.estimator = "kde-adaptive";
else
    info.estimator = "kde";
end
end

% ---------------------------------------------------------------- helpers

function I = mutualInformationBits(X, Y, chunk, adaptive)
n = numel(X);
hx = std(X) / n^(1/6);
hy = std(Y) / n^(1/6);
if hx == 0 || hy == 0
    I = 0;
    return
end

pX = marginalDensity(X, hx, chunk, []);
pY = marginalDensity(Y, hy, chunk, []);
pXY = jointDensity(X, Y, hx, hy, chunk, []);

if adaptive
    % The fixed-bandwidth estimates above are the pilots.
    pX = marginalDensity(X, hx, chunk, abramson(pX));
    pY = marginalDensity(Y, hy, chunk, abramson(pY));
    pXY = jointDensity(X, Y, hx, hy, chunk, abramson(pXY));
end

ok = pXY > 0 & pX > 0 & pY > 0;
I = sum(log2(pXY(ok) ./ (pX(ok) .* pY(ok)))) / n;
end

function p = marginalDensity(v, h, chunk, lambda)
%MARGINALDENSITY Gaussian KDE at the sample points. normpdf inlined.
% LAMBDA, if given, scales the bandwidth of each kernel (one per column).
n = numel(v);
p = zeros(n, 1);
c = 1 / (n * h * sqrt(2*pi));
for lo = 1:chunk:n
    hi = min(lo + chunk - 1, n);
    z = (v(lo:hi) - v.') / h;                 % implicit expansion, no repmat
    if isempty(lambda)
        p(lo:hi) = c * sum(exp(-0.5 * z.^2), 2);
    else
        p(lo:hi) = c * sum(exp(-0.5 * (z ./ lambda.').^2) ./ lambda.', 2);
    end
end
end

function p = jointDensity(X, Y, hx, hy, chunk, lambda)
%JOINTDENSITY Bivariate Gaussian KDE, covariance [hx^2 r*hx*hy; r*hx*hy hy^2].
% Closed form rather than mvnpdf, to avoid the Statistics Toolbox. LAMBDA,
% if given, scales kernel i's covariance by lambda_i^2.
n = numel(X);
r = pearson(X, Y);
r = min(max(r, -0.999999), 0.999999);         % keep the kernel non-singular
omr2 = 1 - r^2;
c = 1 / (n * 2*pi * hx * hy * sqrt(omr2));

p = zeros(n, 1);
for lo = 1:chunk:n
    hi = min(lo + chunk - 1, n);
    u = (X(lo:hi) - X.') / hx;
    v = (Y(lo:hi) - Y.') / hy;
    q = (u.^2 - 2*r*(u.*v) + v.^2) / omr2;
    if isempty(lambda)
        p(lo:hi) = c * sum(exp(-0.5 * q), 2);
    else
        L2 = (lambda.').^2;
        p(lo:hi) = c * sum(exp(-0.5 * q ./ L2) ./ L2, 2);
    end
end
end

function lambda = abramson(pilot)
%ABRAMSON Square-root law: lambda_i = sqrt(g / pilot_i), g the geometric mean.
pilot = max(pilot, realmin);
g = exp(mean(log(pilot)));
lambda = sqrt(g ./ pilot);
end

function r = pearson(a, b)
a = a - mean(a);
b = b - mean(b);
r = (a.' * b) / (norm(a) * norm(b));
end
