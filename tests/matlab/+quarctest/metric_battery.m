function [v, status] = metric_battery(x, fs, M, E, opts)
%METRIC_BATTERY Run the nonlinear metric battery over one series.
%
%   [v, status] = quarctest.metric_battery(x, fs, M, E)
%
%   x   scalar observable, column vector
%   fs  sampling frequency (1 for maps)
%   M   metric policy from quarctest.metric_policy, which decides which
%       columns are computed for this system
%   E   embedding policy from quarctest.embed_policy, which decides what
%       delay and dimension the reconstruction uses
%
%   Returns v, a 1-by-numel(M) row of values aligned with M, and status, a
%   string row recording "ok", "excluded", or the failure for each.
%
%   ESTIMATED AND USED ARE DIFFERENT COLUMNS. ami_delay and fnn_dim are always
%   computed and reported: they are metrics in their own right. What the
%   reconstruction actually uses is embed_delay and embed_dim, which come from
%   E. On flows the two coincide, because a flow's embedding genuinely has to
%   be estimated. On maps they do not: the delay is 1 by construction and the
%   dimension comes from the known state dimension, because AMI and FNN answer
%   a question that a map does not pose. See quarctest.embed_policy.
%
%   Where the embedding is estimated, it is estimated per realization rather
%   than once per system, so the scatter it induces is part of the reported
%   spread -- an analyst re-derives the embedding for each recording, and the
%   variability that matters includes that choice.
%
%   A FAILURE IS A RESULT. Any metric can throw or return a non-finite value
%   on a particular realization -- an FNN search that finds no minimum, a
%   recurrence radius search that fails to converge, a divergence curve with
%   no scaling region. Each is caught, recorded by name in status, and left as
%   NaN in v. The alternative, letting one realization abort the ensemble,
%   would bias the reported distribution towards the systems and initial
%   conditions that happen to be easy.
%
%   See also QUARCTEST.METRIC_POLICY, QUARCTEST.CHARACTERIZE.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    x  (:,1) double
    fs (1,1) double
    M  (1,:) struct
    E  (1,1) struct
    opts.RecTarget (1,1) double = 2.5   % percent recurrence for RQA
    opts.MsScales  (1,1) double = 10    % coarse-graining scales for MSE
    opts.EntDim    (1,1) double = 2     % m for the entropy family
    opts.EntRadius (1,1) double = 0.2   % r for the entropy family, in SD
    opts.PermuDim  (1,1) double = 5     % order for permutation entropy
    opts.MaxLag    (1,1) double = 100   % lag cap for AMI
end

ids = [M.id];
v = nan(1, numel(M));
status = repmat("excluded", 1, numel(M));

want = @(id) any(ids == id & [M.applies]);
put  = @(vv, ss, id, val, st) deal(setAt(vv, ids == id, val), setAt(ss, ids == id, st));

% ---- AMI and FNN, reported as metrics regardless of what drives the battery

amiDelay = NaN; fnnDim = NaN;

try
    amiDelay = ami(x, min(opts.MaxLag, max(2, floor(numel(x)/10))));
    [v, status] = put(v, status, "ami_delay", amiDelay, "ok");
catch err
    [v, status] = put(v, status, "ami_delay", NaN, "fail:" + err.identifier);
end

% FNN needs a delay to search at. It uses the same rule the battery will use,
% so that the reported dimension is the one an analyst following this protocol
% would actually have obtained.
fnnSearchDelay = 1;
if E.delayRule == "ami" && isfinite(amiDelay) && amiDelay >= 1
    fnnSearchDelay = amiDelay;
elseif E.delayRule == "fixed"
    fnnSearchDelay = E.delayValue;
end

try
    fnnDim = fnn(x, fnnSearchDelay, 10, 15, 2, 1);
    [v, status] = put(v, status, "fnn_dim", fnnDim, "ok");
catch err
    [v, status] = put(v, status, "fnn_dim", NaN, "fail:" + err.identifier);
end

% ---- the embedding the rest of the battery actually uses

switch E.delayRule
    case "fixed"
        delay = E.delayValue;
    otherwise
        delay = amiDelay;
end
if ~isfinite(delay) || delay < 1
    delay = 1;      % documented fallback: continue on a stated embedding
end                  % rather than dropping every downstream metric

switch E.dimRule
    case "fixed"
        dim = E.dimValue;
    otherwise
        dim = fnnDim;
        if ~isfinite(dim) || dim < E.stateDim
            dim = E.stateDim + 1;
        end
end

[v, status] = put(v, status, "embed_delay", delay, "ok");
[v, status] = put(v, status, "embed_dim",   dim,   "ok");

% ---- Lyapunov exponents, both methods, in nats per unit time

if want("lyap_wolf")
    [v, status] = tryOne(v, status, ids, "lyap_wolf", ...
        @() lyapunov(x, fs, algorithm="wolf", delay=delay, dim=dim));
end
% Rosenstein is called for its diagnostics as well as its slope: the fitted
% window and its R^2 are what make the exponent interpretable, so they are
% captured from EXTRA in the same call rather than being recomputed or, as
% before, discarded.
if want("lyap_ros") || want("lyap_ros_fit_start") || want("lyap_ros_fit_len") || ...
        want("lyap_ros_fit_r2") || want("lyap_ros_fit_curv") || ...
        want("lyap_ros_fit_maxdev") || want("lyap_ros_fit_runsz")
    fitIds = ["lyap_ros", "lyap_ros_fit_start", "lyap_ros_fit_len", "lyap_ros_fit_r2", ...
              "lyap_ros_fit_curv", "lyap_ros_fit_maxdev", "lyap_ros_fit_runsz"];
    try
        [lam, ex] = lyapunov(x, fs, algorithm="rosenstein", delay=delay, dim=dim);
        lam = firstOf(lam);
        if want("lyap_ros")
            if isfinite(lam)
                [v, status] = put(v, status, "lyap_ros", lam, "ok");
            else
                [v, status] = put(v, status, "lyap_ros", NaN, "nonfinite");
            end
        end
        idx = [];
        if isfield(ex, 'scalingRegion'), idx = ex.scalingRegion(:); end
        if want("lyap_ros_fit_start")
            if isempty(idx)
                [v, status] = put(v, status, "lyap_ros_fit_start", NaN, "nofit");
            else
                [v, status] = put(v, status, "lyap_ros_fit_start", idx(1), "ok");
            end
        end
        if want("lyap_ros_fit_len")
            if isempty(idx)
                [v, status] = put(v, status, "lyap_ros_fit_len", NaN, "nofit");
            else
                [v, status] = put(v, status, "lyap_ros_fit_len", numel(idx), "ok");
            end
        end
        % The divergence curve is fitted against sample index, so the
        % diagnostics are computed on exactly the points lyapunov used.
        ll = struct('r2', NaN, 'curvature', NaN, 'maxDev', NaN, 'runsZ', NaN);
        if ~isempty(idx) && isfield(ex, 'divergence')
            dv = ex.divergence(:);
            ll = quarctest.fit_linearity(idx, dv(idx));
        end
        if want("lyap_ros_fit_r2")
            [v, status] = put(v, status, "lyap_ros_fit_r2", ll.r2, ...
                              ternaryStr(isfinite(ll.r2), "ok", "nofit"));
        end
        if want("lyap_ros_fit_curv")
            [v, status] = put(v, status, "lyap_ros_fit_curv", ll.curvature, ...
                              ternaryStr(isfinite(ll.curvature), "ok", "nofit"));
        end
        if want("lyap_ros_fit_maxdev")
            [v, status] = put(v, status, "lyap_ros_fit_maxdev", ll.maxDev, ...
                              ternaryStr(isfinite(ll.maxDev), "ok", "nofit"));
        end
        if want("lyap_ros_fit_runsz")
            [v, status] = put(v, status, "lyap_ros_fit_runsz", ll.runsZ, ...
                              ternaryStr(isfinite(ll.runsZ), "ok", "nofit"));
        end
    catch err
        for k = 1:numel(fitIds)
            if want(fitIds(k))
                [v, status] = put(v, status, fitIds(k), NaN, "fail:" + err.identifier);
            end
        end
    end
end

% ---- correlation dimension

% As with Rosenstein, corr_dim is called for its diagnostics as well as its
% slope. Its scaling region is chosen by height rather than by linearity, so
% the R^2, the number of bins and the ln(epsilon) span are the only evidence
% that the fitted stretch was straight at all.
cdIds = ["corr_dim", "corr_dim_fit_r2", "corr_dim_fit_len", "corr_dim_fit_span", ...
         "corr_dim_fit_curv", "corr_dim_fit_maxdev", "corr_dim_fit_runsz"];
if any(arrayfun(want, cdIds))
    try
        [cd_, cex] = corr_dim(x, delay, dim, false);
        if want("corr_dim")
            [v, status] = put(v, status, "corr_dim", firstOf(cd_), "ok");
        end
        cl = quarctest.fit_linearity(cex.logEps(cex.idx), cex.logC(cex.idx));
        if want("corr_dim_fit_r2")
            [v, status] = put(v, status, "corr_dim_fit_r2", cl.r2, "ok");
        end
        if want("corr_dim_fit_curv")
            [v, status] = put(v, status, "corr_dim_fit_curv", cl.curvature, "ok");
        end
        if want("corr_dim_fit_maxdev")
            [v, status] = put(v, status, "corr_dim_fit_maxdev", cl.maxDev, "ok");
        end
        if want("corr_dim_fit_runsz")
            [v, status] = put(v, status, "corr_dim_fit_runsz", cl.runsZ, "ok");
        end
        if want("corr_dim_fit_len")
            [v, status] = put(v, status, "corr_dim_fit_len", numel(cex.idx), "ok");
        end
        if want("corr_dim_fit_span")
            le = cex.logEps(:); ii = cex.idx(:);
            if numel(ii) >= 2
                [v, status] = put(v, status, "corr_dim_fit_span", ...
                                  le(ii(end)) - le(ii(1)), "ok");
            else
                [v, status] = put(v, status, "corr_dim_fit_span", NaN, "nofit");
            end
        end
    catch err
        for k = 1:numel(cdIds)
            if want(cdIds(k))
                [v, status] = put(v, status, cdIds(k), NaN, "fail:" + err.identifier);
            end
        end
    end
end

% ---- recurrence quantification, one call feeding every RQA column

rqaIds = ["rqa_radius","rqa_det","rqa_lam","rqa_meanL","rqa_maxL", ...
          "rqa_entL","rqa_entV","rqa_entW"];
if any(ismember(rqaIds, ids(logical([M.applies]))))
    try
        [~, r] = rqa(x, delay, dim, "rec", opts.RecTarget);
        pairs = {"rqa_radius", r.RADIUS; "rqa_det", r.DET; "rqa_lam", r.LAM; ...
                 "rqa_meanL", r.MeanL; "rqa_maxL", r.MaxL; ...
                 "rqa_entL", r.EntrL; "rqa_entV", r.EntrV; "rqa_entW", r.EntrW};
        for k = 1:size(pairs,1)
            if want(pairs{k,1})
                [v, status] = put(v, status, pairs{k,1}, double(pairs{k,2}), "ok");
            end
        end
    catch err
        for k = 1:numel(rqaIds)
            if want(rqaIds(k))
                [v, status] = put(v, status, rqaIds(k), NaN, "fail:" + err.identifier);
            end
        end
    end
end

% ---- entropies

if want("ent_samp")
    [v, status] = tryOne(v, status, ids, "ent_samp", ...
        @() ent_samp(x, opts.EntDim, opts.EntRadius));
end
if want("ent_ap")
    [v, status] = tryOne(v, status, ids, "ent_ap", ...
        @() ent_ap(x, opts.EntDim, opts.EntRadius));
end
if want("ent_permu")
    [v, status] = tryOne(v, status, ids, "ent_permu", ...
        @() ent_permu(x, opts.PermuDim, 1));
end

if want("ent_ms_s1") || want("ent_ms_ci")
    try
        rcmse = ent_ms_plus(x, opts.MsScales, opts.EntDim, opts.EntRadius);
        rcmse = rcmse(:);
        if want("ent_ms_s1")
            [v, status] = put(v, status, "ent_ms_s1", rcmse(1), "ok");
        end
        if want("ent_ms_ci")
            % Complexity index: the summed area under the MSE curve, the
            % standard single-number summary of a multiscale profile.
            [v, status] = put(v, status, "ent_ms_ci", sum(rcmse(isfinite(rcmse))), "ok");
        end
    catch err
        if want("ent_ms_s1")
            [v, status] = put(v, status, "ent_ms_s1", NaN, "fail:" + err.identifier);
        end
        if want("ent_ms_ci")
            [v, status] = put(v, status, "ent_ms_ci", NaN, "fail:" + err.identifier);
        end
    end
end

% ---- detrended fluctuation analysis

if want("dfa_alpha")
    [v, status] = tryOne(v, status, ids, "dfa_alpha", @() localDfaAlpha(x));
end
end

% ---------------------------------------------------------------- internals

function [v, status] = tryOne(v, status, ids, id, fn)
try
    val = fn();
    val = firstOf(val);
    if ~isfinite(val)
        v = setAt(v, ids == id, NaN);
        status = setAt(status, ids == id, "nonfinite");
    else
        v = setAt(v, ids == id, val);
        status = setAt(status, ids == id, "ok");
    end
catch err
    v = setAt(v, ids == id, NaN);
    status = setAt(status, ids == id, "fail:" + err.identifier);
end
end

function a = setAt(a, mask, val)
a(mask) = val;
end

function s = ternaryStr(cond, a, b)
if cond, s = a; else, s = b; end
end

function y = firstOf(z)
z = z(:);
if isempty(z), y = NaN; else, y = double(z(1)); end
end

function alpha = localDfaAlpha(x)
%LOCALDFAALPHA DFA over a decade of scales, returning the fitted exponent.
n = numel(x);
scales = unique(round(logspace(log10(16), log10(max(32, floor(n/8))), 12)));
[~, ~, alpha] = dfa(x, scales, 1, false);
alpha = firstOf(alpha);
end
