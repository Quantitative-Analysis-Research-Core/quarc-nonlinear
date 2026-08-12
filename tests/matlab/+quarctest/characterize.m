function [S, per] = characterize(opts)
%CHARACTERIZE Metric battery over an ensemble of realizations, per system.
%
%   S = quarctest.characterize()
%   S = quarctest.characterize(R=1000, N=4000)
%   S = quarctest.characterize(Fast=true)
%   [S, per] = quarctest.characterize(Systems=["lorenz","henon"], R=100)
%
%   S    summary table, one row per system and metric, with centre, spread,
%        and the published reference where one exists
%   per  per-realization table, one row per system, realization and metric
%
%   WHAT THIS IS. For each system in the catalogue this draws R initial
%   conditions i.i.d. around the published one (Sprott's own x0 is always
%   realization 1), generates a series from each, runs the metric battery, and
%   reports the distribution of every metric across the ensemble. It is an
%   instrument for characterizing how the library's estimators behave across
%   a diverse set of known systems, and its output is a report.
%
%   IT ASSERTS NOTHING AND CANNOT FAIL. There is no pass criterion anywhere in
%   this file, and it is deliberately not named test*.m, so run_tests does not
%   collect it. Two metrics -- the largest Lyapunov exponent and the
%   correlation dimension -- have published reference values, and for those
%   the summary carries the reference, its published uncertainty where Sprott
%   states one, and the observed-to-reference ratio. Those columns invite a
%   judgement; they do not encode one. For every other metric no reference
%   exists, and inventing a tolerance for them would be a number with no
%   provenance dressed up as a criterion.
%
%   WHAT THE SPREAD MEANS. Realizations differ only in initial condition, so
%   for a dissipative system they are finite windows of the same attractor
%   started at different phases, and the spread is the estimator's sampling
%   variability on that attractor. For a conservative system there is no
%   attractor and the realizations are distinct orbits, so the spread there is
%   a different quantity; the category column carries the distinction into the
%   report. Median and MAD are reported alongside mean and SD because several
%   of these distributions are visibly skewed and a mean is not a summary of
%   them.
%
%   SAMPLING IS HELD FIXED WITHIN A SYSTEM. For flows the decimation factor is
%   derived once, from the published initial condition, and applied to every
%   realization, so that fs is a property of the system rather than of the
%   realization. Without this a perturbed x0 could shift the period estimate
%   by one integrator step and change fs, and the resulting scatter in every
%   rate-dependent metric would be an artefact of the protocol.
%
%   Options
%     R           realizations per system. Default 100.
%     N           samples per realization. Default 4000.
%     Systems     names to run, or "all". Default "all".
%     Fast        true runs a twelve-system subset, two per category.
%     Seed        base seed for the initial-condition draw. Default 20260810.
%     Spread      IC perturbation, as a fraction of max(abs(x0),1). Default 0.01.
%     Parallel    use parfor across realizations. Default true. Without the
%                 Parallel Computing Toolbox parfor runs serially, so this is
%                 safe on base MATLAB; it is exposed mainly to force serial
%                 execution when debugging, since errors inside parfor lose
%                 their stack.
%     Verbose     per-system progress to stdout. Default true.
%     Checkpoint  path to a .mat updated after every system, so a long run
%                 that dies partway can be salvaged. Default "" (off).
%     Metrics     restrict the battery to these metric ids. Default "all".
%                 The embedding is always computed, because everything
%                 downstream depends on it, so a restricted pass reproduces
%                 the same delay and dimension as a full one and its rows can
%                 be merged with a full run made at the same seed. Use it to
%                 add a column without repeating the expensive RQA and
%                 entropy work.
%
%   COST. The battery is roughly 14 s per realization at N=4000 and 4-6 s at
%   N=2000, dominated by the RQA radius search, which rebuilds the recurrence
%   matrix once per bisection step. The full catalogue at R=1000, N=4000 is
%   therefore on the order of 230 core-hours: an overnight job on 32 workers,
%   not something to run casually. Start with Fast=true and R=10.
%
%   See also QUARCTEST.ENSEMBLE_ICS, QUARCTEST.METRIC_POLICY,
%   QUARCTEST.METRIC_BATTERY, QUARCTEST.WRITE_CHARACTERIZATION.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    opts.R        (1,1) double {mustBePositive, mustBeInteger} = 100
    opts.N        (1,1) double {mustBePositive, mustBeInteger} = 4000
    opts.Systems  (1,:) string = "all"
    opts.Fast     (1,1) logical = false
    opts.Seed     (1,1) double = 20260810
    opts.Spread   (1,1) double = 0.01
    opts.Parallel (1,1) logical = true
    opts.Verbose  (1,1) logical = true
    opts.Checkpoint (1,1) string = ""
    opts.Metrics  (1,:) string = "all"
end

c = quarctest.sprott_catalog();
c = c([c.usable]);

if opts.Fast
    c = localSubset(c);
elseif ~(isscalar(opts.Systems) && opts.Systems == "all")
    keep = ismember([c.name], opts.Systems);
    missing = setdiff(opts.Systems, [c.name]);
    if ~isempty(missing)
        error('quarctest:characterize:unknownSystem', ...
              'not in the usable catalogue: %s', strjoin(missing, ', '));
    end
    c = c(keep);
end

rows = {};
perBlocks = {};
t0 = tic;

for i = 1:numel(c)
    sys = c(i);
    M = quarctest.metric_policy(sys);
    E = quarctest.embed_policy(sys);
    if ~(isscalar(opts.Metrics) && opts.Metrics == "all")
        keepM = ismember([M.id], opts.Metrics);
        unknown = setdiff(opts.Metrics, [M.id]);
        if ~isempty(unknown)
            error('quarctest:characterize:unknownMetric', ...
                  'not in the battery: %s', strjoin(unknown, ', '));
        end
        for q = find(~keepM)
            M(q).applies = false;
            M(q).reason  = "not requested in this pass";
        end
    end
    ids = [M.id];
    nM = numel(M);

    % Sampling rate fixed once per system, from the published x0.
    [~, gi0] = quarctest.sprott_series(sys, min(opts.N, 2048));
    decim = gi0.decim;
    if sys.kind == "map", decim = 1; end

    ics = quarctest.ensemble_ics(sys, opts.R, Seed=opts.Seed, Spread=opts.Spread);

    V = nan(opts.R, nM);
    ok = false(opts.R, 1);
    why = strings(opts.R, 1);
    fsUsed = nan(opts.R, 1);

    N = opts.N;
    if opts.Parallel
        parfor r = 1:opts.R
            [V(r,:), ok(r), why(r), fsUsed(r)] = ...
                localOne(sys, ics(:,r), N, decim, M, E); %#ok<PFBNS>
        end
    else
        for r = 1:opts.R
            [V(r,:), ok(r), why(r), fsUsed(r)] = ...
                localOne(sys, ics(:,r), N, decim, M, E);
        end
    end

    if opts.Verbose
        fprintf('%-24s %-20s R=%d  usable=%d  fs=%.4g  %.1fs\n', ...
                sys.name, sys.category, opts.R, sum(ok), ...
                median(fsUsed, 'omitnan'), toc(t0));
    end

    % ---- per-realization rows, built as one block per system
    %
    % Assembled with repmat rather than appended row by row: a full run is
    % 55 systems x R realizations x ~20 metrics, which at R=1000 is over a
    % million rows, and growing a cell array one element at a time reallocates
    % on every step.
    applied = find([M.applies]);
    if ~isempty(applied)
        nA = numel(applied);
        realIdx = repmat((1:opts.R)', nA, 1);
        metIdx  = repelem(ids(applied)', opts.R, 1);
        vals    = reshape(V(:, applied), [], 1);
        okCol   = repmat(ok, nA, 1);
        nRow    = numel(realIdx);
        % The four system columns and the metric name repeat for every row, so
        % they are stored categorical rather than as string arrays: at R=1000
        % across the catalogue that is over a million rows, and five string
        % arrays of that length cost gigabytes to hold identical text.
        perBlocks{end+1} = table( ...
            repmat(categorical(sys.name), nRow, 1), ...
            repmat(categorical(sys.section), nRow, 1), ...
            repmat(categorical(sys.category), nRow, 1), ...
            repmat(categorical(string(sys.kind)), nRow, 1), ...
            realIdx, categorical(metIdx), vals, okCol, ...
            'VariableNames', {'system','section','category','kind', ...
                              'realization','metric','value','seriesUsable'}); %#ok<AGROW>
    end

    % ---- summary rows
    for k = 1:nM
        if ~M(k).applies
            rows{end+1} = localRow(sys, M(k), [], opts, NaN, NaN); %#ok<AGROW>
            continue
        end
        v = V(ok, k);
        v = v(isfinite(v));
        rows{end+1} = localRow(sys, M(k), v, opts, ...
                               median(fsUsed, 'omitnan'), sum(ok)); %#ok<AGROW>
    end

    % ---- checkpoint after every system
    %
    % A full run is a long job. Losing a system's worth of work to a crash is
    % tolerable; losing fifty is not, so the accumulated summary is written
    % after each system and the run can be salvaged from the last one.
    if opts.Checkpoint ~= ""
        ckpt = struct('rows', {rows}, 'systemsDone', i, ...
                      'systemsTotal', numel(c), 'lastSystem', sys.name, ...
                      'R', opts.R, 'N', opts.N, 'seed', opts.Seed, ...
                      'elapsed', toc(t0)); %#ok<NASGU>
        save(opts.Checkpoint, '-struct', 'ckpt');
    end
end

S = localToTable(rows);
if isempty(perBlocks)
    per = table();
else
    per = vertcat(perBlocks{:});
end

if opts.Verbose
    fprintf('\ncharacterize: %d systems, R=%d, N=%d, %.1f s total\n', ...
            numel(c), opts.R, opts.N, toc(t0));
end
end

% ---------------------------------------------------------------- internals

function [v, ok, why, fs] = localOne(sys, ic, N, decim, M, E)
%LOCALONE One realization: series from this IC, then the battery.
v = nan(1, numel(M));
ok = false; why = ""; fs = NaN;
s = sys;
s.x0 = ic(:);
try
    if s.kind == "map"
        [x, gi] = quarctest.sprott_series(s, N);
    else
        [x, gi] = quarctest.sprott_series(s, N, Decim=decim);
    end
catch err
    why = "series:" + err.identifier;
    return
end
if gi.degenerate
    why = "degenerate";
    return
end
fs = gi.fs;
v = quarctest.metric_battery(x, fs, M, E);
ok = true;
end

function row = localRow(sys, m, vals, opts, fs, nUsable)
%LOCALROW One summary row: centre, spread, and the reference if there is one.
if isempty(vals)
    n = 0; med = NaN; mad_ = NaN; mu = NaN; sd = NaN; lo = NaN; hi = NaN;
else
    n = numel(vals);
    med = median(vals);
    mad_ = median(abs(vals - med));     % MAD about the median, unscaled
    mu = mean(vals);
    sd = std(vals);
    lo = min(vals);
    hi = max(vals);
end

refVal = NaN; refErr = NaN; refSrc = "";
switch m.reference
    case "lambda"
        refVal = sys.lambda;
        refErr = NaN;                    % tier says exact or numerical
        refSrc = "Sprott A, lambda (" + sys.tier + ")";
    case "d2"
        refVal = sys.d2;
        refErr = sys.d2_err;
        refSrc = "Sprott A, D2";
end

if isfinite(refVal) && isfinite(med) && refVal ~= 0
    ratio = med / refVal;
else
    ratio = NaN;
end

row = {sys.name, sys.section, sys.category, string(sys.kind), ...
       m.id, m.label, m.units, ...
       string(ternary(m.applies, "characterized", "excluded")), ...
       m.reason, opts.N, opts.R, nUsable, fs, n, ...
       med, mad_, mu, sd, lo, hi, ...
       refVal, refErr, refSrc, ratio, ...
       string(ternary(m.reference == "", "none", "published")), ...
       opts.Seed, opts.Spread};
end

function T = localToTable(rows)
if isempty(rows), T = table(); return, end
A = vertcat(rows{:});
T = cell2table(A, 'VariableNames', { ...
    'system','section','category','kind', ...
    'metric','label','units','role','excludedBecause', ...
    'N','R','usableRealizations','fs','n', ...
    'median','mad','mean','sd','min','max', ...
    'reference','referenceUncertainty','referenceSource','medianOverReference', ...
    'referenceStatus','seed','spread'});
end

function c = localSubset(c)
%LOCALSUBSET Twelve systems, two per category, for a fast pass.
% Two per category, chosen to span the appendix rather than to be easy:
% conservative systems and driven flows are the hardest cases in the set and
% a fast subset that omitted them would flatter the library.
want = ["logistic","cusp", ...                  % noninvertible maps
        "henon","ikeda", ...                    % dissipative maps
        "chirikov","arnold_cat", ...            % conservative maps
        "driven_vdp","ueda", ...                % driven flows
        "lorenz","rossler", ...                 % autonomous flows
        "nose_hoover","henon_heiles"];          % conservative flows
keep = ismember([c.name], want);
if ~any(keep)
    keep = false(1, numel(c));
    keep(round(linspace(1, numel(c), min(10, numel(c))))) = true;
end
c = c(keep);
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
