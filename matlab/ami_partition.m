function [tau, curve, info] = ami_partition(x, maxlag, opts)
%AMI_PARTITION Average mutual information versus lag, adaptive partition.
%   TAU = AMI_PARTITION(X,maxlag) returns the lag of the first local minimum
%   of the average mutual information of X with itself over lags 0:maxlag,
%   estimated on an adaptive partition of the plane that is refined only
%   where the data have structure (Fraser and Swinney 1986; Darbellay and
%   Vajda 1999).
%
%   [TAU,CURVE] = AMI_PARTITION(X,maxlag) also returns the AMI curve as an
%   (maxlag+1)-by-2 array of [lag, ami], in bits.
%
%   [TAU,CURVE,INFO] = AMI_PARTITION(X,maxlag) also returns a struct with
%   fields allMinima, usedFallback, fractionLag, samplesPerLag, alpha and
%   estimator.
%
%   ___ = AMI_PARTITION(X,maxlag,alpha=A) sets the significance level of the
%   chi-square test that decides whether a cell is split. Default 0.05.
%   Smaller values split less and give a coarser, lower estimate.
%
%   ___ = AMI_PARTITION(X,maxlag,fraction=F) sets the fallback threshold, as
%   a fraction of AMI at lag 0, used when the curve has no local minimum.
%   Default 0.2.
%
%   Input Arguments
%      X       time series, real column vector, no NaN
%      maxlag  maximum lag, positive integer, less than numel(X)
%
%   Notes
%   Both coordinates are replaced by their ranks, so every marginal is
%   uniform and the estimate does not depend on the units or the shape of
%   the distribution of X. Starting from the whole plane, each cell is cut
%   at the median rank of each coordinate into four. The cut is kept if a
%   chi-square test (3 degrees of freedom) rejects a uniform spread over the
%   four sub-cells, and the test is repeated on each sub-cell. The root is
%   always cut. A cell that is not cut contributes
%
%      (n_c/N) * log2(n_c*N / (w_x*w_y))
%
%   where n_c is its count and w_x, w_y its widths in ranks, which are the
%   marginal counts because ranks are a permutation (Darbellay and Vajda
%   1999). Ties are broken by position.
%
%   Every lag uses the same number of pairs, N-maxlag, as in AMI_HISTOGRAM.
%   At lag 0 the points lie on the diagonal and the partition refines until
%   cells hold about four points, so AMI(0) is roughly log2(N)-2 bits. That
%   value is a property of the estimator, as it is for every estimator of a
%   continuous variable; it matters only through the Fraction fallback.
%
%   Runs on base MATLAB. No toolbox required.
%
%   Examples
%      tau = ami_partition(x, 50);
%      [tau, curve] = ami_partition(x, 50, alpha=0.01);
%
%   References
%      Fraser, A. M. and Swinney, H. L. (1986). Independent coordinates for
%      strange attractors from mutual information. Physical Review A, 33(2),
%      1134-1140.
%
%      Darbellay, G. A. and Vajda, I. (1999). Estimation of the information
%      by an adaptive partitioning of the observation space. IEEE
%      Transactions on Information Theory, 45(4), 1315-1321.
%
%   See also AMI, AMI_HISTOGRAM, AMI_KDE, AMI_KNN, AMI_COPULA.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    x  (:,1) double {mustBeNonempty}
    maxlag  (1,1) double {mustBePositive, mustBeInteger}
    opts.alpha    (1,1) double {mustBeGreaterThan(opts.alpha, 0), mustBeLessThan(opts.alpha, 1)} = 0.05
    opts.fraction (1,1) double {mustBePositive} = 0.2
end

if anynan(x)
    error('ami_partition:nanInput', ...
          'x contains NaN. Mutual information is undefined for missing data; remove or impute first.');
end
N = numel(x);
if maxlag >= N
    error('ami_partition:lagTooLarge', ...
          'maxlag (%d) must be smaller than the number of samples (%d).', maxlag, N);
end

critical = chi2Critical3(opts.alpha);
m = N - maxlag;

curve = zeros(maxlag + 1, 2);
for lag = 0:maxlag
    a = x(1:m);
    b = x((1 + lag):(m + lag));
    curve(lag + 1, :) = [lag, partitionBits(a, b, critical)];
end

[tau, info] = ami_first_minimum(curve, opts.fraction, "ami_partition");
info.samplesPerLag = m;
info.alpha = opts.alpha;
info.estimator = "adaptive";
end

% ---------------------------------------------------------------- helpers

function I = partitionBits(a, b, critical)
%PARTITIONBITS Darbellay-Vajda estimate of I(A;B) in bits.
n = numel(a);
rx = rankOf(a);
ry = rankOf(b);

% Explicit stack rather than recursion: depth grows with log2(N), but an
% explicit stack keeps the cost of a cell to the cell.
bounds = [1 n 1 n];               % [xlo xhi ylo yhi], inclusive, in ranks
members = {(1:n).'};
isRoot = true;
I = 0;

while ~isempty(bounds)
    c = bounds(end, :);
    idx = members{end};
    bounds(end, :) = [];
    members(end) = [];

    nc = numel(idx);
    if nc == 0
        continue
    end
    wx = c(2) - c(1) + 1;
    wy = c(4) - c(3) + 1;

    split = false;
    if wx >= 2 && wy >= 2
        xm = c(1) + floor(wx/2) - 1;          % last rank of the lower half
        ym = c(3) + floor(wy/2) - 1;
        lx = rx(idx) <= xm;
        ly = ry(idx) <= ym;
        q = [sum(lx & ly); sum(lx & ~ly); sum(~lx & ly); sum(~lx & ~ly)];

        % Expected counts under a uniform spread, by area: halves are unequal
        % when a width is odd.
        fx = (xm - c(1) + 1) / wx;
        fy = (ym - c(3) + 1) / wy;
        E = nc * [fx*fy; fx*(1-fy); (1-fx)*fy; (1-fx)*(1-fy)];
        split = isRoot || sum((q - E).^2 ./ E) > critical;
    end
    isRoot = false;

    if split
        bounds = [bounds;                                   %#ok<AGROW>
                  c(1)   xm   c(3)   ym;
                  c(1)   xm   ym+1   c(4);
                  xm+1   c(2) c(3)   ym;
                  xm+1   c(2) ym+1   c(4)];
        members = [members; {idx(lx & ly); idx(lx & ~ly); ...  %#ok<AGROW>
                             idx(~lx & ly); idx(~lx & ~ly)}];
    else
        I = I + (nc / n) * log2(nc * n / (wx * wy));
    end
end
end

function r = rankOf(v)
%RANKOF Ranks 1..n, ties broken by position. sort is stable.
[~, order] = sort(v);
r = zeros(numel(v), 1);
r(order) = 1:numel(v);
end

function c = chi2Critical3(alpha)
%CHI2CRITICAL3 Upper-alpha quantile of chi-square with 3 degrees of freedom.
% Closed-form CDF, so chi2inv (Statistics Toolbox) is not needed:
%   F(x) = erf(sqrt(x/2)) - sqrt(2x/pi) * exp(-x/2)
F = @(x) erf(sqrt(x/2)) - sqrt(2*x/pi) .* exp(-x/2);
c = fzero(@(x) F(x) - (1 - alpha), [1e-9, 200]);
end
