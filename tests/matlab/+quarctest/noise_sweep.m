function [S, per] = noise_sweep(opts)
%NOISE_SWEEP Estimator accuracy against observational noise, with and without re-embedding.
%
%   [S, per] = quarctest.noise_sweep()
%   [S, per] = quarctest.noise_sweep(Systems=["Lorenz","Rossler"], R=10)
%   [S, per] = quarctest.noise_sweep(Noise=[0 0.01 0.05], Checkpoint=c)
%
%   THE QUESTION. How much observational noise does each exponent estimator
%   survive, and how much of the damage is the estimator's rather than the
%   embedding's? Noise inflates the AMI first minimum and the FNN dimension, so
%   an analyst following this protocol on a noisy series embeds it differently
%   than they would the clean one. Those are two different failures and they
%   are separated here by running every noisy series twice:
%
%     arm "fixed"    delay and dimension taken from the CLEAN series and held.
%                    What the estimator alone does to a noisier signal.
%     arm "reembed"  delay and dimension re-estimated on the noisy series, the
%                    way a practitioner would. What a real analysis produces.
%
%   Both arms see the SAME noise draw on the SAME realization, so a difference
%   between them is the embedding and nothing else. Neither arm is the honest
%   answer on its own: "reembed" is what a user gets, "fixed" is what the
%   estimator is responsible for, and the gap between them is the result.
%
%   The embedding parameters are carried as outcomes too (ami_delay, fnn_dim,
%   embed_delay, embed_dim), because "noise inflates the delay" is a claim this
%   sweep can measure rather than assume.
%
%   WHAT WOULD DISTINGUISH THE EXPLANATIONS. If the two arms track each other,
%   the embedding is not what breaks and the estimators are simply noise
%   sensitive. If "fixed" holds up while "reembed" degrades, the estimators are
%   sound and the protocol's embedding step is the weak link -- which would be
%   an argument for embedding on a smoothed copy, not for a different estimator.
%
%   NOISE MODEL. Additive Gaussian, scaled to each series' own standard
%   deviation, so a level means the same thing on systems with different units.
%   It is measurement noise on the stored observable: the trajectory is not
%   perturbed, only the view of it. Noise injected into the integration is a
%   different experiment and needs its own export.
%
%   Level 0 is the clean series. Both arms are identical there by construction
%   -- there is nothing to re-embed differently -- and it is computed once and
%   recorded under both arms so each arm carries a complete curve. That cell is
%   also the check on this harness: at R=100 it must reproduce the lyap_wolf,
%   lyap_ros and corr_dim rows of tests/reports/characterization_dysts.csv.
%
%   Options
%     Manifest    dysts manifest to read (default: the committed one)
%     Systems     "all" or a list of system names
%     R           realizations per system (default 25)
%     Noise       noise levels as a fraction of series SD
%                 (default [0 0.005 0.01 0.02 0.05 0.10])
%     Seed        base seed for the noise draws (default 20260810)
%     Parallel    parfor across realizations (default true)
%     Verbose     one progress line per system (default true)
%     Checkpoint  .mat path, rewritten after every system
%
%   Returns
%     S    one row per (system, noise, arm, metric): n, median, mad,
%          reference and ratio
%     per  one row per (system, realization, noise, arm, metric)
%
%   See also QUARCTEST.CHARACTERIZE_SERIES, QUARCTEST.EVOLVE_SWEEP.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    opts.Manifest   (1,1) string  = ""
    opts.Systems    (1,:) string  = "all"
    opts.R          (1,1) double {mustBePositive, mustBeInteger} = 25
    opts.Noise      (1,:) double {mustBeNonnegative} = [0 0.005 0.01 0.02 0.05 0.10]
    opts.Seed       (1,1) double  = 20260810
    opts.Parallel   (1,1) logical = true
    opts.Verbose    (1,1) logical = true
    opts.Checkpoint (1,1) string  = ""
    % Which metric family to sweep. "lyapunov" is the committed 16 September
    % record; "rqa" is the pass that record deferred until line_hist counted
    % lines correctly. The embedding metrics come along with either.
    opts.Family     (1,1) string {mustBeMember(opts.Family, ["lyapunov", "rqa", "all"])} = "lyapunov"
end

quarctest.require_library();
c = quarctest.dysts_catalog(opts.Manifest);
if ~(isscalar(opts.Systems) && opts.Systems == "all")
    keep = ismember([c.name], opts.Systems);
    missing = setdiff(opts.Systems, [c.name]);
    if ~isempty(missing)
        error('quarctest:noiseSweep:unknownSystem', ...
              'not in the manifest: %s', strjoin(missing, ', '));
    end
    c = c(keep);
end

% The metrics this sweep reports. The embedding pair is nearly free -- the
% battery computes ami and fnn whatever drives it -- and they are the outcome
% of the second experiment, so they are kept rather than discarded.
% The first sweep (16 September) covered the Lyapunov estimators and the
% correlation dimension and left the RQA family out on purpose: line_hist's
% counting fix was not yet on main, and sweeping before it landed would have
% measured the bug rather than the metric. That fix is on main now, so the RQA
% family is swept here. The families run as separate passes rather than one,
% because the Lyapunov pass is the committed record and would be identical.
EMBED = ["ami_delay", "fnn_dim", "embed_delay", "embed_dim"];
switch opts.Family
    case "lyapunov"
        WANT = [EMBED, "lyap_wolf", "lyap_ros", "corr_dim"];
    case "rqa"
        WANT = [EMBED, "rqa_radius", "rqa_det", "rqa_lam", "rqa_maxL", ...
                "rqa_maxL_tw", "rqa_meanL", "rqa_entL", "rqa_entV", "rqa_entW"];
    case "all"
        WANT = [EMBED, "lyap_wolf", "lyap_ros", "corr_dim", ...
                "rqa_radius", "rqa_det", "rqa_lam", "rqa_maxL", ...
                "rqa_maxL_tw", "rqa_meanL", "rqa_entL", "rqa_entV", "rqa_entW"];
end

levels = opts.Noise(:)';
nL = numel(levels);
perBlocks = {};
t0 = tic;

for i = 1:numel(c)
    sys = c(i);
    M = quarctest.metric_policy(sys);
    E = quarctest.embed_policy(sys);

    % Same subsetting idiom characterize_series uses: the battery still returns
    % every id, the ones not asked for are marked so they are not computed.
    keepM = ismember([M.id], WANT);
    for q = find(~keepM)
        M(q).applies = false;
        M(q).reason  = "not requested in this pass";
    end
    ids = [M.id];
    wantIdx = find(ismember(ids, WANT));

    X = localRead(sys);
    R = min(opts.R, size(X, 1));
    X = X(1:R, :);
    fs = sys.fs;

    % One cell per realization, each holding a (level x arm x metric) page.
    cells = cell(R, 1);
    base = opts.Seed + i * 100000;

    if opts.Parallel
        parfor r = 1:R
            cells{r} = localOne(X(r,:), fs, M, E, levels, base + r); %#ok<PFBNS>
        end
    else
        for r = 1:R
            cells{r} = localOne(X(r,:), fs, M, E, levels, base + r);
        end
    end

    perBlocks{end+1} = localPer(sys, cells, ids, wantIdx, levels); %#ok<AGROW>

    if opts.Verbose
        usable = sum(cellfun(@(z) ~isempty(z) && z.ok, cells));
        fprintf('%-28s R=%d usable=%d levels=%d  %.1f s\n', ...
                sys.name, R, usable, nL, toc(t0));
    end

    if opts.Checkpoint ~= ""
        ckpt = struct('perBlocks', {perBlocks}, 'systemsDone', i, ...
                      'systemsTotal', numel(c), 'lastSystem', sys.name, ...
                      'elapsed', toc(t0)); %#ok<NASGU>
        save(opts.Checkpoint, '-struct', 'ckpt');
    end
end

if isempty(perBlocks)
    per = table(); S = table(); return
end
per = vertcat(perBlocks{:});
S = localSummarize(per, c);

if opts.Verbose
    fprintf('\nnoise_sweep: %d systems, R=%d, %d levels, %.1f s\n', ...
            numel(c), opts.R, nL, toc(t0));
end
end

% ---------------------------------------------------------------- internals

function X = localRead(sys)
%LOCALREAD The ensemble as written by export_dysts.py. Identical to the reader
%   in characterize_series: numpy writes row-major, so read [N, R] and transpose.
fid = fopen(sys.file, 'r');
if fid < 0
    error('quarctest:noiseSweep:noFile', 'cannot open %s', sys.file);
end
cl = onCleanup(@() fclose(fid));
raw = fread(fid, [sys.N, sys.R], sys.dtype);
X = double(raw');
end

function out = localOne(x, fs, M, E, levels, seed)
%LOCALONE One realization: every noise level, both arms.
out = struct('ok', false, 'V', []);
x = x(:);
if any(~isfinite(x)) || std(x) <= 0
    return
end

nL = numel(levels);
nM = numel(M);
V = nan(nL, 2, nM);            % level x arm(1=fixed,2=reembed) x metric
sd = std(x);

% The clean pass fixes the embedding the "fixed" arm will carry. embed_delay
% and embed_dim are reported by the battery whatever the rules say, so this
% reads back the protocol's own choice rather than recomputing it here.
clean = quarctest.metric_battery(x, fs, M, E);
ids = [M.id];
delay0 = clean(ids == "embed_delay");
dim0   = clean(ids == "embed_dim");

Efix = E;
Efix.delayRule = "fixed"; Efix.delayValue = delay0;
Efix.dimRule   = "fixed"; Efix.dimValue   = dim0;

% One stream per realization, drawn in level order, so a run is reproducible
% and both arms of a level see the same numbers.
rs = RandStream('threefry', 'Seed', mod(seed, 2^31));

for k = 1:nL
    if levels(k) == 0
        % Nothing to perturb and nothing to re-embed: one computation, recorded
        % under both arms so each carries a complete curve from 0 upward.
        V(k,1,:) = clean;
        V(k,2,:) = clean;
        continue
    end
    y = x + (levels(k) * sd) * randn(rs, numel(x), 1);
    V(k,1,:) = quarctest.metric_battery(y, fs, M, Efix);
    V(k,2,:) = quarctest.metric_battery(y, fs, M, E);
end

out.ok = true;
out.V  = V;
end

function T = localPer(sys, cells, ids, wantIdx, levels)
%LOCALPER Long-form rows for one system.
arms = ["fixed", "reembed"];
R = numel(cells);
rows = {};
for r = 1:R
    z = cells{r};
    if isempty(z) || ~z.ok, continue, end
    for k = 1:numel(levels)
        for a = 1:2
            for q = wantIdx
                rows{end+1} = {string(sys.name), string(sys.category), ...
                               r, levels(k), arms(a), ids(q), ...
                               z.V(k, a, q)}; %#ok<AGROW>
            end
        end
    end
end
if isempty(rows), T = table(); return, end
T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'system','category','realization','noise','arm','metric','value'});
T.system = categorical(T.system);
T.category = categorical(T.category);
T.arm = categorical(T.arm);
T.metric = categorical(T.metric);
end

function S = localSummarize(per, c)
%LOCALSUMMARIZE Median and MAD per (system, noise, arm, metric), with the
%   published reference where one exists.
lam = containers.Map(cellstr(string([c.name])), num2cell([c.lambda]));
d2  = containers.Map(cellstr(string([c.name])), num2cell([c.d2]));

[g, sysv, noisev, armv, metv] = findgroups(per.system, per.noise, per.arm, per.metric);
med = splitapply(@(v) median(v(isfinite(v))), per.value, g);
mad_ = splitapply(@(v) median(abs(v(isfinite(v)) - median(v(isfinite(v))))), per.value, g);
n   = splitapply(@(v) sum(isfinite(v)), per.value, g);

S = table(sysv, noisev, armv, metv, n, med, mad_, ...
    'VariableNames', {'system','noise','arm','metric','n','median','mad'});

S.reference = nan(height(S), 1);
for j = 1:height(S)
    name = char(string(S.system(j)));
    switch string(S.metric(j))
        case {"lyap_wolf", "lyap_ros"}
            S.reference(j) = lam(name);
        case "corr_dim"
            S.reference(j) = d2(name);
    end
end
S.medianOverReference = S.median ./ S.reference;
S.medianOverReference(~isfinite(S.reference) | S.reference == 0) = NaN;
S = sortrows(S, {'system','metric','arm','noise'});
end
