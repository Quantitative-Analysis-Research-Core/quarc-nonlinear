function [tau, curve, info] = ami_copula(x, maxlag, opts)
%AMI_COPULA Average mutual information versus lag, Gaussian copula estimator.
%   TAU = AMI_COPULA(X,maxlag) returns the lag of the first local minimum of
%   the Gaussian-copula mutual information of X with itself over lags
%   0:maxlag (Ince et al. 2017).
%
%   [TAU,CURVE] = AMI_COPULA(X,maxlag) also returns the AMI curve as an
%   (maxlag+1)-by-2 array of [lag, ami], in bits.
%
%   [TAU,CURVE,INFO] = AMI_COPULA(X,maxlag) also returns a struct with
%   fields allMinima, usedFallback, fractionLag, samplesPerLag and estimator.
%
%   ___ = AMI_COPULA(X,maxlag,fraction=F) sets the fallback threshold, as a
%   fraction of AMI at lag 0, used when the curve has no local minimum.
%   Default 0.2.
%
%   Input Arguments
%      X       time series, real column vector, no NaN
%      maxlag  maximum lag, positive integer, less than numel(X)
%
%   Notes
%   Each coordinate is replaced by its ranks and mapped through the inverse
%   standard normal CDF, which makes both marginals exactly Gaussian. The
%   mutual information of a bivariate Gaussian with correlation r is
%
%      I = -0.5 * log2(1 - r^2)
%
%   and that is the estimate. It is exact when the dependence between X(t)
%   and X(t+lag) has a Gaussian copula, and a lower bound otherwise. It is
%   fast, O(maxlag*N*log N), and has very low variance.
%
%   It sees only dependence that a correlation of normal scores can see.
%   Where X(t+lag) depends on X(t) non-monotonically -- a sine at a quarter
%   period, where the pairs trace a circle -- r is near zero and so is the
%   estimate, although the true AMI is large. Its first minimum therefore
%   sits near the first zero of the rank autocorrelation, which is a
%   different criterion from Fraser and Swinney's. Compare it with another
%   estimator before relying on it for a delay.
%
%   At lag 0, r = 1 and the true value is infinite; r is capped at
%   1 - 1e-6, as in AMI_KDE, giving about 9.6 bits.
%
%   Every lag uses the same number of pairs, N-maxlag, as in AMI_HISTOGRAM.
%   Ties are broken by position.
%
%   Runs on base MATLAB. No toolbox required.
%
%   Examples
%      tau = ami_copula(x, 50);
%
%   References
%      Ince, R. A. A., Giordano, B. L., Kayser, C., Rousselet, G. A.,
%      Gross, J. and Schyns, P. G. (2017). A statistical framework for
%      neuroimaging data analysis based on mutual information estimated via
%      a Gaussian copula. Human Brain Mapping, 38(3), 1541-1573.
%
%   See also AMI, AMI_HISTOGRAM, AMI_KDE, AMI_KNN, AMI_PARTITION.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    x  (:,1) double {mustBeNonempty}
    maxlag  (1,1) double {mustBePositive, mustBeInteger}
    opts.fraction (1,1) double {mustBePositive} = 0.2
end

if anynan(x)
    error('ami_copula:nanInput', ...
          'x contains NaN. Mutual information is undefined for missing data; remove or impute first.');
end
N = numel(x);
if maxlag >= N
    error('ami_copula:lagTooLarge', ...
          'maxlag (%d) must be smaller than the number of samples (%d).', maxlag, N);
end
m = N - maxlag;

curve = zeros(maxlag + 1, 2);
for lag = 0:maxlag
    a = normalScores(x(1:m));
    b = normalScores(x((1 + lag):(m + lag)));
    r = (a.' * b) / (norm(a) * norm(b));
    r = min(max(r, -0.999999), 0.999999);
    curve(lag + 1, :) = [lag, -0.5 * log2(1 - r^2)];
end

[tau, info] = ami_first_minimum(curve, opts.fraction, "ami_copula");
info.samplesPerLag = m;
info.estimator = "copula";
end

% ---------------------------------------------------------------- helpers

function g = normalScores(v)
%NORMALSCORES Ranks mapped to standard normal quantiles. norminv inlined
% through erfcinv, so the Statistics Toolbox is not needed.
n = numel(v);
[~, order] = sort(v);
u = zeros(n, 1);
u(order) = (1:n).' / (n + 1);
g = -sqrt(2) * erfcinv(2 * u);
end
