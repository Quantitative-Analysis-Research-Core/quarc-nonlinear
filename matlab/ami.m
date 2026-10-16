function [tau, curve, info] = ami(x, maxlag, opts)
%AMI Average mutual information versus lag.
%   TAU = AMI(X,maxlag) returns the lag of the first local minimum of the average
%   mutual information of X with itself, computed over lags 0:maxlag. TAU is the
%   embedding delay for phase space reconstruction.
%
%   [TAU,CURVE] = AMI(X,maxlag) also returns the AMI curve as an (maxlag+1)-by-2 array
%   whose columns are [lag, ami]. AMI is in bits.
%
%   [TAU,CURVE,INFO] = AMI(X,maxlag) also returns a struct with fields:
%      allMinima     all local minima found, as [lag, ami] rows
%      usedFallback  true if no local minimum existed and the fraction
%                    threshold was used instead
%      fractionLag   lag at which AMI first fell below Fraction*AMI(0)
%      bins          number of bins used ("histogram" only)
%      k             number of neighbours ("knn" only)
%      estimator     estimator that ran
%      algorithm     value of Algorithm as supplied
%
%   ___ = AMI(X,maxlag,algorithm=ALG) selects the estimator:
%      "histogram"  (default) equal-width joint histogram. Fast.
%                   Aliases: "stergiou".
%      "kde"        bivariate Gaussian kernel density estimate. Slower,
%                   O(maxlag*N^2), but less biased. Aliases: "thomas".
%      "kde-adaptive"  as "kde", with each kernel's width set by the local
%                   density (Abramson's square-root law). About 2x "kde".
%      "adaptive"   adaptive partition of the rank plane, refined where a
%                   chi-square test finds structure (Fraser and Swinney;
%                   Darbellay and Vajda). Aliases: "partition".
%      "knn"        Kraskov-Stoegbauer-Grassberger nearest-neighbour
%                   estimator. No bins or bandwidth; near-unbiased.
%                   O(maxlag*N^2). Aliases: "ksg".
%      "copula"     Gaussian copula. Fastest and lowest variance, but blind
%                   to non-monotonic dependence: see AMI_COPULA before using
%                   it to choose a delay. Aliases: "gcmi".
%
%   ___ = AMI(X,maxlag,bins=B) sets the number of histogram bins per axis.
%   Default selects B by Scott's rule. Used only by "histogram".
%
%   ___ = AMI(X,maxlag,k=K) sets the number of neighbours. Default 10. Used
%   only by "knn".
%
%   ___ = AMI(X,maxlag,fraction=F) sets the fallback threshold, as a fraction of
%   AMI at lag 0, used when the curve has no local minimum. Default 0.2.
%
%   Input Arguments
%      X  time series, real column vector, no NaN
%      maxlag  maximum lag, positive integer, less than numel(X)
%
%   Notes
%   The estimators disagree in absolute value; the histogram is biased
%   upward and the bias grows with bin count. For a bivariate Gaussian with
%   correlation r the true AMI is -0.5*log2(1-r^2); on an AR(1) series with
%   phi = 0.7 at lag 1 the exact value is 0.4857 bits, against 0.5619
%   (histogram) and 0.5411 (kde). Report which estimator produced a
%   published tau. Options for the individual estimators that AMI does not
%   pass through (alpha, chunk) are set by calling them directly.
%
%   Examples
%      % Embedding delay for a Lorenz series
%      tau = ami(x, 50);
%
%      % Same series, kernel density estimator
%      tau = ami(x, 50, algorithm="kde");
%
%      % Fixed bin count, and inspect the curve
%      [tau, curve] = ami(x, 50, bins=64);
%      plot(curve(:,1), curve(:,2))
%
%   References
%      Fraser, A. M. and Swinney, H. maxlag. (1986). Independent coordinates for
%      strange attractors from mutual information. Physical Review A, 33(2),
%      1134-1140.
%
%      Thomas, R. D., Moses, N. C., Semple, E. A. and Strang, A. J. (2014).
%      An efficient algorithm for the computation of average mutual
%      information. Behavior Research Methods.
%
%      Abarbanel, H. D. I., Brown, R., Sidorowich, J. J. and Tsimring, maxlag. S.
%      (1993). The analysis of observed chaotic data in physical systems.
%      Reviews of Modern Physics, 65(4), 1331-1392.
%
%   See also AMI_HISTOGRAM, AMI_KDE, AMI_PARTITION, AMI_KNN, AMI_COPULA,
%   LYE_R, LYE_W.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    x  (:,1) double {mustBeNonempty}
    maxlag  (1,1) double {mustBePositive, mustBeInteger}
    opts.algorithm (1,1) string = "histogram"
    opts.bins      double {mustBeScalarOrEmpty, mustBePositive, mustBeInteger} = []
    opts.k         double {mustBeScalarOrEmpty, mustBePositive, mustBeInteger} = []
    opts.fraction  (1,1) double {mustBePositive} = 0.2
end

algorithm = lower(opts.algorithm);
if ~isempty(opts.bins) && ~ismember(algorithm, ["histogram", "stergiou"])
    warning('ami:binsIgnored', ...
            'bins applies to algorithm="histogram" and is ignored for "%s".', algorithm);
end
if ~isempty(opts.k) && ~ismember(algorithm, ["knn", "ksg"])
    warning('ami:kIgnored', ...
            'k applies to algorithm="knn" and is ignored for "%s".', algorithm);
end

switch algorithm

    case {"histogram", "stergiou"}
        args = {"Fraction", opts.fraction};
        if ~isempty(opts.bins)
            args = [args, {"Bins", opts.bins}];
        end
        [tau, curve, info] = ami_histogram(x, maxlag, args{:});

    case {"kde", "thomas"}
        [tau, curve, info] = ami_kde(x, maxlag, "Fraction", opts.fraction);

    case "kde-adaptive"
        [tau, curve, info] = ami_kde(x, maxlag, "Fraction", opts.fraction, ...
                                     "bandwidth", "adaptive");

    case {"adaptive", "partition"}
        [tau, curve, info] = ami_partition(x, maxlag, "Fraction", opts.fraction);

    case {"knn", "ksg"}
        args = {"Fraction", opts.fraction};
        if ~isempty(opts.k)
            args = [args, {"k", opts.k}];
        end
        [tau, curve, info] = ami_knn(x, maxlag, args{:});

    case {"copula", "gcmi"}
        [tau, curve, info] = ami_copula(x, maxlag, "Fraction", opts.fraction);

    otherwise
        error('ami:unknownAlgorithm', ...
              ['Unknown algorithm "%s". Supported: "histogram" (alias "stergiou"), ' ...
               '"kde" (alias "thomas"), "kde-adaptive", "adaptive" (alias ' ...
               '"partition"), "knn" (alias "ksg"), "copula" (alias "gcmi").'], ...
              opts.algorithm);
end

info.algorithm = algorithm;
end
