function [tau, curve, info] = ami_knn(x, maxlag, opts)
%AMI_KNN Average mutual information versus lag, k-nearest-neighbour estimator.
%   TAU = AMI_KNN(X,maxlag) returns the lag of the first local minimum of the
%   average mutual information of X with itself over lags 0:maxlag, using
%   the Kraskov-Stoegbauer-Grassberger estimator (algorithm 1).
%
%   [TAU,CURVE] = AMI_KNN(X,maxlag) also returns the AMI curve as an
%   (maxlag+1)-by-2 array of [lag, ami], in bits.
%
%   [TAU,CURVE,INFO] = AMI_KNN(X,maxlag) also returns a struct with fields
%   allMinima, usedFallback, fractionLag, samplesPerLag, k and estimator.
%
%   ___ = AMI_KNN(X,maxlag,k=K) sets the number of neighbours. Default 10.
%   Larger K lowers variance and, in principle, raises bias. Kraskov et al.
%   recommend 2-4 for estimating a value; the default here is larger
%   because AMI is used to find a minimum, and that needs a smooth curve.
%   On a period-40 sine with 15% noise (N = 2000, ten draws) the first
%   minimum fell at lags 4-9 with K = 3 and 8-11 with K = 10, the expected
%   lag being 10, while the AR(1) error at lag 1 was unchanged
%   (0.4945 vs 0.4931 bits, exact 0.4857).
%
%   ___ = AMI_KNN(X,maxlag,fraction=F) sets the fallback threshold, as a
%   fraction of AMI at lag 0, used when the curve has no local minimum.
%   Default 0.2.
%
%   ___ = AMI_KNN(X,maxlag,chunk=C) sets how many rows of the pairwise
%   distance matrix are evaluated at once. Default 512. It changes peak
%   memory, not the result.
%
%   Input Arguments
%      X       time series, real column vector, no NaN
%      maxlag  maximum lag, positive integer, less than numel(X)
%
%   Notes
%   For each pair (x_i, y_i), eps_i is the max-norm distance to its k-th
%   nearest neighbour in the joint space, and n_x(i), n_y(i) count the
%   points strictly within eps_i in each marginal. Then, in nats,
%
%      I = psi(k) + psi(N) - < psi(n_x + 1) + psi(n_y + 1) >
%
%   converted here to bits. The estimate has no bins and no bandwidth and
%   is close to unbiased for independent data, so unlike AMI_HISTOGRAM it
%   can return small negative values where the true AMI is zero.
%
%   The series is standardised, then perturbed by 1e-10 standard deviations
%   to break ties, as Kraskov et al. advise. The perturbation is not random:
%   sample i moves by 1e-10*(mod(i*g, 1) - 0.5), g the golden-ratio
%   conjugate, an evenly spread sequence that every language computes
%   bit for bit. Repeated calls return identical results, the global random
%   stream is untouched, and ports reproduce the result exactly.
%
%   Every lag uses the same number of pairs, N-maxlag, as in AMI_HISTOGRAM.
%   Cost is O(maxlag*N^2) by brute-force neighbour search; peak memory is
%   O(Chunk*N).
%
%   Runs on base MATLAB. No toolbox required.
%
%   Examples
%      tau = ami_knn(x, 50);
%      [tau, curve] = ami_knn(x, 50, k=4);   % closer to the paper
%
%   References
%      Kraskov, A., Stoegbauer, H. and Grassberger, P. (2004). Estimating
%      mutual information. Physical Review E, 69(6), 066138.
%
%   See also AMI, AMI_HISTOGRAM, AMI_KDE, AMI_PARTITION, AMI_COPULA.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    x  (:,1) double {mustBeNonempty}
    maxlag  (1,1) double {mustBePositive, mustBeInteger}
    opts.k        (1,1) double {mustBePositive, mustBeInteger} = 10
    opts.fraction (1,1) double {mustBePositive} = 0.2
    opts.chunk    (1,1) double {mustBePositive, mustBeInteger} = 512
end

if anynan(x)
    error('ami_knn:nanInput', ...
          'x contains NaN. Mutual information is undefined for missing data; remove or impute first.');
end
N = numel(x);
if maxlag >= N
    error('ami_knn:lagTooLarge', ...
          'maxlag (%d) must be smaller than the number of samples (%d).', maxlag, N);
end
m = N - maxlag;
if opts.k >= m
    error('ami_knn:kTooLarge', ...
          'k (%d) must be smaller than the number of pairs per lag (%d).', opts.k, m);
end

s = std(x);
if s > 0
    x = (x - mean(x)) / s;
end
x = x + 1e-10 * (mod((1:N).' * 0.6180339887498949, 1) - 0.5);

curve = zeros(maxlag + 1, 2);
for lag = 0:maxlag
    a = x(1:m);
    b = x((1 + lag):(m + lag));
    curve(lag + 1, :) = [lag, ksgBits(a, b, opts.k, opts.chunk)];
end

[tau, info] = ami_first_minimum(curve, opts.fraction, "ami_knn");
info.samplesPerLag = m;
info.k = opts.k;
info.estimator = "knn";
end

% ---------------------------------------------------------------- helpers

function I = ksgBits(a, b, k, chunk)
%KSGBITS Kraskov-Stoegbauer-Grassberger algorithm 1, in bits.
n = numel(a);
nx = zeros(n, 1);
ny = zeros(n, 1);
for lo = 1:chunk:n
    hi = min(lo + chunk - 1, n);
    rows = (lo:hi).';
    dx = abs(a(rows) - a.');                  % implicit expansion
    dy = abs(b(rows) - b.');
    d = max(dx, dy);
    d(sub2ind(size(d), (1:numel(rows)).', rows)) = Inf;   % exclude self
    e = mink(d, k, 2);
    e = e(:, end);
    % Strictly closer than eps_i; the self term (distance 0) is counted by
    % the comparison and removed.
    nx(rows) = sum(dx < e, 2) - 1;
    ny(rows) = sum(dy < e, 2) - 1;
end
I = (psi(k) + psi(n) - mean(psi(nx + 1) + psi(ny + 1))) / log(2);
end
