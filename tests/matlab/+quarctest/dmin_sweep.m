function [S, per] = dmin_sweep(opts)
%DMIN_SWEEP The DET family against its minimum line length.
%
%   [S, per] = quarctest.dmin_sweep()
%   [S, per] = quarctest.dmin_sweep(R=100, Dmin=[2 5 10 20 40])
%
%   For every dysts system, every realization: reconstruct the phase space
%   exactly as the battery does, find the radius for the protocol recurrence
%   target exactly as rqa's "rec" branch does, and take the diagonal line
%   histogram ONCE. dmin is a post-hoc threshold on that histogram, so the
%   whole curve -- DET, MeanL, MaxL, EntrL at every requested dmin -- costs
%   one O(N^2) pass instead of one per value. The dmin=2 column doubles as
%   the corrected protocol value, superseding the report rows computed
%   before line_hist counted lines exactly.
%
%   THE QUESTION. rqa's dmin defaults to 2, and the export samples 40 points
%   per dominant period, so a "determinism-qualifying" diagonal spans 1/20th
%   of an orbit -- tangential motion alone hands that to nearly every
%   recurrent point, and DET sits within a fraction of a percent of 100
%   across the whole catalogue. The recurrence-target sweep showed the
%   radius is the wrong knob: loosening it saturates DET harder. This sweep
%   asks whether dmin is the right one; at dmin=40 a line must survive one
%   full orbital period.
%
%   Options
%     Manifest    path to the dysts manifest.json. Default tests/reports/dysts.
%     Systems     names to run, or "all".
%     R           realizations per system, capped at what was exported.
%     Dmin        minimum line lengths, in samples. Default [2 5 10 20 40].
%     RecTarget   percent recurrence for the radius search. Default 2.5.
%     MaxLag      lag cap for AMI, as in metric_battery. Default 100.
%     Parallel    parfor across realizations. Default true.
%     Verbose     one line per system. Default true.
%     Checkpoint  .mat rewritten after every system.
%
%   Returns
%     S     one row per system and dmin: n, median/mad of DET, MeanL, MaxL,
%           EntrL across the ensemble, and the median radius.
%     per   one row per realization and dmin, for distribution work.
%
%   See also QUARCTEST.CHARACTERIZE_SERIES, QUARCTEST.EVOLVE_SWEEP, RQA.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    opts.Manifest   (1,1) string = ""
    opts.Systems    (1,:) string = "all"
    opts.R          (1,1) double {mustBePositive, mustBeInteger} = 100
    opts.Dmin       (1,:) double {mustBePositive, mustBeInteger} = [2 5 10 20 40]
    opts.RecTarget  (1,1) double = 2.5
    opts.MaxLag     (1,1) double = 100
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
        error('quarctest:dminSweep:unknownSystem', ...
              'not in the manifest: %s', strjoin(missing, ', '));
    end
    c = c(keep);
end

dmins = opts.Dmin;
nD = numel(dmins);
perBlocks = {};
t0 = tic;

for i = 1:numel(c)
    sys = c(i);
    E = quarctest.embed_policy(sys);
    ts = tic;

    X = localRead(sys);
    R = min(opts.R, size(X, 1));
    X = X(1:R, :);

    DET = nan(R, nD); MEANL = nan(R, nD); MAXL = nan(R, nD); ENTL = nan(R, nD);
    RAD = nan(R, 1); ok = false(R, 1);
    rt = opts.RecTarget; ml = opts.MaxLag;

    if opts.Parallel
        parfor r = 1:R
            [DET(r,:), MEANL(r,:), MAXL(r,:), ENTL(r,:), RAD(r), ok(r)] = ...
                localOne(X(r,:)', E, dmins, rt, ml); %#ok<PFBNS>
        end
    else
        for r = 1:R
            [DET(r,:), MEANL(r,:), MAXL(r,:), ENTL(r,:), RAD(r), ok(r)] = ...
                localOne(X(r,:)', E, dmins, rt, ml);
        end
    end

    if opts.Verbose
        fprintf('%-28s R=%d  usable=%d  %.1fs\n', sys.name, R, sum(ok), toc(ts));
    end

    rIdx = repmat((1:R)', nD, 1);
    dCol = repelem(dmins(:), R, 1);
    perBlocks{end+1} = table( ...
        repmat(categorical(sys.name), R*nD, 1), ...
        repmat(categorical(sys.category), R*nD, 1), ...
        rIdx, dCol, ...
        reshape(DET, [], 1), reshape(MEANL, [], 1), ...
        reshape(MAXL, [], 1), reshape(ENTL, [], 1), ...
        repmat(RAD, nD, 1), repmat(ok, nD, 1), ...
        'VariableNames', {'system','category','realization','dmin', ...
                          'det','meanL','maxL','entL','radius','seriesUsable'}); %#ok<AGROW>

    if opts.Checkpoint ~= ""
        ckpt = struct('perBlocks', {perBlocks}, 'lastSystem', sys.name, ...
                      'systemsDone', i, 'systemsTotal', numel(c), ...
                      'elapsed', toc(t0));
        save(opts.Checkpoint, '-struct', 'ckpt');
    end
end

per = vertcat(perBlocks{:});
S = localSummarize(per, dmins);
fprintf('dmin_sweep: %d systems, %d dmin values, %.1f s\n', ...
        numel(c), nD, toc(t0));
end

% ---------------------------------------------------------------- internals

function [det_, meanL, maxL, entL, rad, ok] = localOne(x, E, dmins, recTarget, maxLag)
% The reconstruction and radius search reproduce metric_battery calling
% rqa(x, fs, algorithm..., "rec", target) step for step, so a dmin=2 column
% here is the protocol value, not an approximation of it.
nD = numel(dmins);
det_ = nan(1, nD); meanL = nan(1, nD); maxL = nan(1, nD); entL = nan(1, nD);
rad = NaN; ok = false;
if any(~isfinite(x)) || std(x) <= 0
    return
end

amiDelay = NaN;
try
    amiDelay = ami(x, min(maxLag, max(2, floor(numel(x)/10))));
catch
end
switch E.delayRule
    case "fixed"
        delay = E.delayValue;
    otherwise
        delay = amiDelay;
end
if ~isfinite(delay) || delay < 1
    delay = 1;
end

fnnDim = NaN;
try
    fnnDim = fnn(x, delay, 10, 15, 2, 1);
catch
end
switch E.dimRule
    case "fixed"
        dim = E.dimValue;
    otherwise
        dim = fnnDim;
        if ~isfinite(dim) || dim < E.stateDim
            dim = E.stateDim + 1;
        end
end

try
    y = psr(x, delay, dim);
    a = -abs(pdist2(y, y));
    [rec, dh, ~, rad] = set_radius(y, a, 0.01, 0.5, recTarget, 'rqa', 20);
    dh = dh(:)'; dh = dh(dh > 0);
    if rec <= 0 || isempty(dh)
        return
    end
    for k = 1:nD
        q = dh(dh >= dmins(k));
        det_(k) = 100 * sum(q) / sum(dh);
        if isempty(q)
            continue          % DET is a true 0; the line stats stay NaN
        end
        meanL(k) = mean(q);
        maxL(k) = max(q);
        count = histcounts(q, (min(q) - 0.5):(max(q) + 0.5));
        p = count(count > 0) / sum(count);
        entL(k) = -sum(p .* log2(p));
    end
    ok = true;
catch
end
end

function X = localRead(sys)
fid = fopen(sys.file, 'r');
if fid < 0
    error('quarctest:dminSweep:noSeries', 'cannot open %s', sys.file);
end
raw = fread(fid, [sys.N, sys.R], sys.dtype);
fclose(fid);
X = raw';
end

function S = localSummarize(per, dmins)
rows = {};
for sysName = unique(per.system, 'stable')'
    % The category comes from the system, not from the usable subset: a
    % system whose every realization failed (PanXuZhou cannot reach the
    % recurrence target) still gets its n=0 rows rather than crashing the
    % summary after the whole sweep has run.
    cat_ = per.category(find(per.system == sysName, 1));
    for d = dmins
        q = per(per.system == sysName & per.dmin == d & per.seriesUsable, :);
        v = q(isfinite(q.det), :);
        rows{end+1} = {sysName, cat_, d, height(v), ...
            median(v.det), mad_(v.det), median(v.meanL), mad_(v.meanL), ...
            median(v.maxL), mad_(v.maxL), median(v.entL), mad_(v.entL), ...
            median(v.radius)}; %#ok<AGROW>
    end
end
S = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'system','category','dmin','n', ...
     'det_median','det_mad','meanL_median','meanL_mad', ...
     'maxL_median','maxL_mad','entL_median','entL_mad','radius_median'});
end

function m = mad_(v)
v = v(isfinite(v));
if isempty(v), m = NaN; else, m = median(abs(v - median(v))); end
end
