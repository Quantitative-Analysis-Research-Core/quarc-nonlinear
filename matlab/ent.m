function [value, extra] = ent(x, opts)
%ENT Entropy of a time series.
%   VALUE = ENT(X) returns the sample entropy of X.
%
%   VALUE is a scalar for every algorithm except "multiscale", which returns
%   a column of entropies, one per coarse-graining scale.
%
%   [VALUE,EXTRA] = ENT(...) also returns a struct with the algorithm that
%   ran, the arguments it used, and any secondary outputs the method
%   produces.
%
%   ___ = ENT(X,algorithm=ALG) selects the method:
%      "sample"      (default) sample entropy. Richman & Moorman (2000).
%      "approximate" approximate entropy. Pincus (1991). Biased by
%                    self-matches; prefer "sample" for new work.
%      "permutation" permutation entropy of ordinal patterns. Bandt &
%                    Pompe (2002). Needs no radius.
%      "symbolic"    normalised corrected Shannon entropy of a symbol
%                    sequence. X must be binary; threshold first, e.g.
%                    x = double(x > median(x)).
%      "multiscale"  refined composite multiscale sample entropy and its
%                    relatives, over coarse-graining scales 1:scale.
%
%   ___ = ENT(X,y=Y,algorithm=ALG) computes a cross entropy between X and Y:
%      "sample"      cross sample entropy
%      "approximate" cross approximate entropy
%
%   ___ = ENT(X,dim=D) sets the template length. Default 2.
%   ___ = ENT(X,radius=R) sets the match tolerance as a proportion of the
%   standard deviation. Default 0.2. Unused by "permutation" and "symbolic".
%   ___ = ENT(X,delay=T) sets the ordinal-pattern delay. "permutation" only,
%   default 1.
%   ___ = ENT(X,scale=S) sets the greatest coarse-graining scale.
%   "multiscale" only, default 20. This is a scale count, not a delay.
%
%   Input Arguments
%      X  time series, real column vector, no NaN
%
%   Notes
%   dim and radius are m and r in Richman & Moorman (2000). Sample entropy
%   counts template matches of length dim and dim+1 over the same N-dim
%   templates, so the ratio is a conditional probability.
%
%   Approximate entropy counts self-matches, which biases it toward
%   regularity and makes it depend on series length. It is kept for
%   comparison with older results; sample entropy is the better estimator.
%
%   Permutation and symbolic entropy work from ordinal patterns and symbol
%   sequences, so they take no radius and are insensitive to a monotone
%   rescaling of the signal.
%
%   Examples
%      value = ent(x);                              % sample entropy
%      value = ent(x, algorithm="permutation", dim=3);
%      value = ent(x, algorithm="multiscale", scale=10);
%      value = ent(x, y=y, algorithm="sample");     % cross sample entropy
%
%      % All five on one series
%      for a = ["sample" "approximate" "permutation" "symbolic" "multiscale"]
%          fprintf('%-12s %.4f\n', a, ent(x, algorithm=a));
%      end
%
%   References
%      Richman, J. S. and Moorman, J. R. (2000). Physiological time-series
%      analysis using approximate entropy and sample entropy. American
%      Journal of Physiology, 278(6), H2039-H2049.
%
%      Pincus, S. M. (1991). Approximate entropy as a measure of system
%      complexity. PNAS, 88(6), 2297-2301.
%
%      Bandt, C. and Pompe, B. (2002). Permutation entropy: a natural
%      complexity measure for time series. Physical Review Letters, 88(17),
%      174102.
%
%   See also ENT_SAMP, ENT_AP, ENT_PERMU, ENT_SYMBOLIC, ENT_MS_PLUS,
%   ENT_XSAMP, ENT_XAP, AMI, LYAPUNOV.

arguments
    x               (:,1) double {mustBeNonempty}
    opts.algorithm  (1,1) string = "sample"
    opts.y          double = []
    opts.dim        (1,1) double {mustBePositive, mustBeInteger} = 2
    opts.radius     (1,1) double {mustBePositive} = 0.2
    opts.delay      (1,1) double {mustBePositive, mustBeInteger} = 1
    opts.scale      (1,1) double {mustBePositive, mustBeInteger} = 20
end

if anynan(x)
    error('ent:nanInput', ...
        'x contains NaN. Entropy is undefined for missing data; remove or impute first.');
end
if numel(x) <= opts.dim + 1
    error('ent:tooShort', ...
        'Series of %d samples is too short for dim = %d.', numel(x), opts.dim);
end

isCross = ~isempty(opts.y);
if isCross
    y = opts.y(:);
    if anynan(y)
        error('ent:nanInput', 'y contains NaN.');
    end
    if numel(y) ~= numel(x)
        error('ent:lengthMismatch', ...
            'x has %d samples and y has %d; a cross entropy needs equal lengths.', ...
            numel(x), numel(y));
    end
end

extra = struct('algorithm', lower(opts.algorithm), 'dim', opts.dim, ...
               'cross', isCross);

switch lower(opts.algorithm)

    case "sample"
        extra.radius = opts.radius;
        if isCross
            value = ent_xsamp(x, y, opts.dim, opts.radius, 1);
            extra.estimator = "cross sample";
        else
            value = ent_samp(x, opts.dim, opts.radius);
            extra.estimator = "sample";
        end

    case "approximate"
        extra.radius = opts.radius;
        if isCross
            value = ent_xap(x, y, opts.dim, opts.radius, 1);
            extra.estimator = "cross approximate";
        else
            value = ent_ap(x, opts.dim, opts.radius);
            extra.estimator = "approximate";
        end

    case "permutation"
        localRejectCross(isCross, "permutation");
        [value, counts] = ent_permu(x, opts.dim, opts.delay);
        extra.delay = opts.delay;
        extra.patternCounts = counts;
        extra.estimator = "permutation";

    case "symbolic"
        localRejectCross(isCross, "symbolic");
        u = unique(x);
        if ~all(ismember(u, [0 1]))
            error('ent:notBinary', ...
                ['algorithm "symbolic" needs a binary series; x has %d ' ...
                 'distinct values in [%g, %g]. Threshold it first, for ' ...
                 'example x = double(x > median(x)).'], ...
                numel(u), min(u), max(u));
        end
        value = ent_symbolic(x, opts.dim);
        extra.estimator = "symbolic";

    case "multiscale"
        localRejectCross(isCross, "multiscale");
        [rcmse, cmse, mse, msfe, gmse] = ent_ms_plus(x, opts.scale, opts.dim, opts.radius);
        value = rcmse;
        extra.radius = opts.radius;
        extra.scale = opts.scale;
        extra.rcmse = rcmse;
        extra.cmse  = cmse;
        extra.mse   = mse;
        extra.msfe  = msfe;
        extra.gmse  = gmse;
        extra.estimator = "refined composite multiscale sample";

    otherwise
        error('ent:unknownAlgorithm', ...
            ['Unknown algorithm "%s". Supported: "sample", "approximate", ' ...
             '"permutation", "symbolic", "multiscale".'], opts.algorithm);
end
end

% ---------------------------------------------------------------- helper

function localRejectCross(isCross, name)
if isCross
    error('ent:noCrossVariant', ...
        ['algorithm "%s" has no cross-entropy form. Supply y only with ' ...
         '"sample" or "approximate".'], name);
end
end
