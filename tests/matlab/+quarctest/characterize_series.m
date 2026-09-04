function [S, per] = characterize_series(opts)
%CHARACTERIZE_SERIES Metric battery over ensembles that were generated elsewhere.
%
%   S = quarctest.characterize_series()
%   [S, per] = quarctest.characterize_series(R=100, Systems=["Lorenz","Chua"])
%
%   Same battery, same policies, same output schema as quarctest.characterize.
%   The difference is where the series come from: characterize integrates each
%   system itself, this reads ensembles written by python/export_dysts.py.
%
%   WHY A SECOND ENTRY POINT. The dysts vector fields are Python, and the
%   published exponents this suite compares against were computed by that same
%   Python. Re-implementing 129 systems in MATLAB would put the data and its
%   reference values on different provenances, and a transcription slip in a
%   coefficient produces a perfectly plausible chaotic series that nothing
%   downstream would flag. Generating the trajectories where the references
%   were generated keeps the two tied together.
%
%   The protocol on the Python side mirrors quarctest.ensemble_ics: realization
%   1 is the published initial condition, the rest are drawn on the attractor
%   and relaxed, and sampling is fixed at 40 points per dominant period, which
%   is what quarctest.sprott_series decimates towards. So an fs here means what
%   an fs there means.
%
%   Options
%     Manifest    path to the dysts manifest. Default tests/reports/dysts.
%     Systems     names to run, or "all".
%     R           realizations per system, capped at what was exported.
%     Metrics     restrict the battery, as in quarctest.characterize.
%     Parallel    parfor across realizations. Default true.
%     Checkpoint  .mat rewritten after every system.
%
%   See also QUARCTEST.CHARACTERIZE, QUARCTEST.DYSTS_CATALOG.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    opts.Manifest   (1,1) string = ""
    opts.Systems    (1,:) string = "all"
    opts.R          (1,1) double {mustBePositive, mustBeInteger} = 100
    opts.Metrics    (1,:) string = "all"
    opts.Parallel   (1,1) logical = true
    opts.Verbose    (1,1) logical = true
    opts.Checkpoint (1,1) string = ""
end

quarctest.require_library();
c = quarctest.dysts_catalog(opts.Manifest);
if ~(isscalar(opts.Systems) && opts.Systems == "all")
    keep = ismember([c.name], opts.Systems);
    missing = setdiff(opts.Systems, [c.name]);
    if ~isempty(missing)
        error('quarctest:characterizeSeries:unknownSystem', ...
              'not in the manifest: %s', strjoin(missing, ', '));
    end
    c = c(keep);
end

rows = {}; perBlocks = {};
t0 = tic;

for i = 1:numel(c)
    sys = c(i);
    M = quarctest.metric_policy(sys);
    E = quarctest.embed_policy(sys);
    if ~(isscalar(opts.Metrics) && opts.Metrics == "all")
        keepM = ismember([M.id], opts.Metrics);
        unknown = setdiff(opts.Metrics, [M.id]);
        if ~isempty(unknown)
            error('quarctest:characterizeSeries:unknownMetric', ...
                  'not in the battery: %s', strjoin(unknown, ', '));
        end
        for q = find(~keepM)
            M(q).applies = false;
            M(q).reason  = "not requested in this pass";
        end
    end
    ids = [M.id];
    nM = numel(M);

    X = localRead(sys);                 % R-by-N, dtype named in the manifest
    R = min(opts.R, size(X, 1));
    X = X(1:R, :);
    fs = sys.fs;

    V = nan(R, nM);
    ok = false(R, 1);

    if opts.Parallel
        parfor r = 1:R
            [V(r,:), ok(r)] = localOne(X(r,:), fs, M, E); %#ok<PFBNS>
        end
    else
        for r = 1:R
            [V(r,:), ok(r)] = localOne(X(r,:), fs, M, E);
        end
    end

    if opts.Verbose
        fprintf('%-28s %-24s R=%d  usable=%d  fs=%.4g  %.1fs\n', ...
                sys.name, sys.category, R, sum(ok), fs, toc(t0));
    end

    applied = find([M.applies]);
    if ~isempty(applied)
        nA = numel(applied);
        realIdx = repmat((1:R)', nA, 1);
        metIdx  = repelem(ids(applied)', R, 1);
        vals    = reshape(V(:, applied), [], 1);
        okCol   = repmat(ok, nA, 1);
        nRow    = numel(realIdx);
        perBlocks{end+1} = table( ...
            repmat(categorical(sys.name), nRow, 1), ...
            repmat(categorical(sys.section), nRow, 1), ...
            repmat(categorical(sys.category), nRow, 1), ...
            repmat(categorical(string(sys.kind)), nRow, 1), ...
            realIdx, categorical(metIdx), vals, okCol, ...
            'VariableNames', {'system','section','category','kind', ...
                              'realization','metric','value','seriesUsable'}); %#ok<AGROW>
    end

    for k = 1:nM
        if ~M(k).applies
            rows{end+1} = localRow(sys, M(k), [], R, sys.N, NaN, NaN); %#ok<AGROW>
            continue
        end
        v = V(ok, k); v = v(isfinite(v));
        rows{end+1} = localRow(sys, M(k), v, R, sys.N, fs, sum(ok)); %#ok<AGROW>
    end

    if opts.Checkpoint ~= ""
        ckpt = struct('rows', {rows}, 'systemsDone', i, ...
                      'systemsTotal', numel(c), 'lastSystem', sys.name, ...
                      'elapsed', toc(t0)); %#ok<NASGU>
        save(opts.Checkpoint, '-struct', 'ckpt');
    end
end

S = localToTable(rows);
if isempty(perBlocks), per = table(); else, per = vertcat(perBlocks{:}); end

if opts.Verbose
    fprintf('\ncharacterize_series: %d systems, R=%d, %.1f s\n', ...
            numel(c), opts.R, toc(t0));
end
end

% ---------------------------------------------------------------- internals

function X = localRead(sys)
fid = fopen(sys.file, 'r');
if fid < 0
    error('quarctest:characterizeSeries:noFile', 'cannot open %s', sys.file);
end
cl = onCleanup(@() fclose(fid));
raw = fread(fid, [sys.N, sys.R], sys.dtype);   % written row-major from numpy
X = double(raw');
end

function [v, ok] = localOne(x, fs, M, E)
v = nan(1, numel(M));
ok = false;
x = x(:);
if any(~isfinite(x)) || std(x) <= 0
    return
end
v = quarctest.metric_battery(x, fs, M, E);
ok = true;
end

function row = localRow(sys, m, vals, R, N, fs, nUsable)
if isempty(vals)
    n = 0; med = NaN; mad_ = NaN; mu = NaN; sd = NaN; lo = NaN; hi = NaN;
else
    n = numel(vals); med = median(vals); mad_ = median(abs(vals - med));
    mu = mean(vals); sd = std(vals); lo = min(vals); hi = max(vals);
end

refVal = NaN; refErr = NaN; refSrc = "";
switch m.reference
    case "lambda"
        refVal = sys.lambda; refSrc = "dysts, tangent-space spectrum";
    case "d2"
        refVal = sys.d2;     refSrc = "dysts, correlation dimension";
end
if isfinite(refVal) && isfinite(med) && refVal ~= 0
    ratio = med/refVal;
else
    ratio = NaN;
end

row = {sys.name, sys.section, sys.category, string(sys.kind), ...
       m.id, m.label, m.units, ...
       string(ternary(m.applies, "characterized", "excluded")), ...
       m.reason, N, R, nUsable, fs, n, ...
       med, mad_, mu, sd, lo, hi, ...
       refVal, refErr, refSrc, ratio, ...
       string(ternary(m.reference == "", "none", "published")), ...
       NaN, NaN};
end

function T = localToTable(rows)
if isempty(rows), T = table(); return, end
T = cell2table(vertcat(rows{:}), 'VariableNames', { ...
    'system','section','category','kind', ...
    'metric','label','units','role','excludedBecause', ...
    'N','R','usableRealizations','fs','n', ...
    'median','mad','mean','sd','min','max', ...
    'reference','referenceUncertainty','referenceSource','medianOverReference', ...
    'referenceStatus','seed','spread'});
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
