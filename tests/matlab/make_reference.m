function make_reference()
%MAKE_REFERENCE Record MATLAB's answers on the shared fixtures.
%
%   matlab -batch "addpath('tests/matlab'); make_reference"
%
%   Writes tests/fixtures/matlab_reference.json, which IS COMMITTED. The
%   R suite (rpkg/tests) compares against this file, so cross-language equivalence is
%   checked in CI without needing a MATLAB licence on the runner.
%
%   Regenerate deliberately, never as part of a test run: if a test could
%   rewrite its own expected values, a port that drifted would simply update
%   the reference and keep passing. Fixtures and reference are inputs to the
%   tests, not outputs of them.
%
%   Only DETERMINISTIC quantities go in here. Anything driven by rand() cannot
%   be compared value-by-value across languages -- MATLAB and NumPy do not
%   share a generator -- and is tested by contract instead (see
%   quarctest.surrogateContract), not by equality.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(fullfile(repo, 'matlab'));
addpath(here);
dbclear all

fx = fullfile(repo, 'tests', 'fixtures');
ref = struct();

series = {'white_4096', 'fgn_H30_4096', 'fgn_H50_4096', 'fgn_H85_4096'};
for i = 1:numel(series)
    y = readmatrix(fullfile(fx, [series{i} '.csv']));
    n = numel(y);
    sc = unique(round(logspace(log10(16), log10(floor(n/8)), 18)));
    [~, fl, alpha] = dfa(y, sc, 1, false);
    ref.dfa.(series{i}) = struct( ...
        'scales', sc, 'fluctuation', fl(:)', 'alpha', alpha);
end

y = readmatrix(fullfile(fx, 'ar1_phi70_512.csv'));
ref.ent_samp.ar1_phi70_512 = ent_samp(y, 2, 0.2);

y = readmatrix(fullfile(fx, 'lorenz_2048.csv'));
ref.ami_thomas.lorenz_2048 = localFirstMinAmi(y);

% ---- References for the R port (rpkg/). Each records the full output, not
% just a summary, so a port that agrees on the summary by accident still fails.

short = {'ar1_phi70_512', 'henon_512', 'pseudoperiodic_300'};
for i = 1:numel(short)
    y = readmatrix(fullfile(fx, [short{i} '.csv']));
    ref.ent_samp_dims.(short{i}) = [ent_samp(y, 1, 0.2), ent_samp(y, 2, 0.2), ...
                                    ent_samp(y, 3, 0.15)];
    ref.ent_ap.(short{i}) = [ent_ap(y, 1, 0.2), ent_ap(y, 2, 0.2), ent_ap(y, 3, 0.15)];
end

y = readmatrix(fullfile(fx, 'lorenz_2048.csv'));
Y = psr(y, 10, 3);
ref.psr.lorenz_2048_d10_m3 = struct('size', size(Y), 'first', Y(1,:), ...
    'last', Y(end,:), 'colsum', sum(Y, 1));
[pe3, h3] = ent_permu(y, 3, 1);
[pe4, h4] = ent_permu(y, 4, 2);
ref.ent_permu.lorenz_2048 = struct('dim3_delay1', pe3, 'hist3', h3, ...
                                   'dim4_delay2', pe4, 'hist4', h4);
y = readmatrix(fullfile(fx, 'white_4096.csv'));
[pe, h] = ent_permu(y, 4, 1);
ref.ent_permu.white_4096 = struct('dim4_delay1', pe, 'hist4', h);

amiSeries = {'lorenz_2048', 'white_4096', 'ar1_phi70_512'};
for i = 1:numel(amiSeries)
    y = readmatrix(fullfile(fx, [amiSeries{i} '.csv']));
    [tau, curve, info] = ami_histogram(y, 20);
    ref.ami_histogram.(amiSeries{i}) = struct('tau', tau, 'curve', curve(:,2)', ...
        'bins', info.bins);
    [tau, curve] = ami_histogram(y, 20, bins=16);
    ref.ami_histogram_bins16.(amiSeries{i}) = struct('tau', tau, 'curve', curve(:,2)');
end
y = readmatrix(fullfile(fx, 'lorenz_2048.csv'));
[tau, curve] = ami_kde(y, 12);
ref.ami_kde.lorenz_2048 = struct('tau', tau, 'curve', curve(:,2)');
y = readmatrix(fullfile(fx, 'ar1_phi70_512.csv'));
[tau, curve] = ami_kde(y, 8);
ref.ami_kde.ar1_phi70_512 = struct('tau', tau, 'curve', curve(:,2)');

% The estimators added after the first R batch. knn breaks ties with a
% deterministic sequence rather than random noise, so it is comparable too.
amiNew = {'lorenz_2048', 'ar1_phi70_512', 'henon_512'};
for i = 1:numel(amiNew)
    y = readmatrix(fullfile(fx, [amiNew{i} '.csv']));
    [tau, curve] = ami_partition(y, 15);
    ref.ami_partition.(amiNew{i}) = struct('tau', tau, 'curve', curve(:,2)');
    [tau, curve] = ami_partition(y, 15, alpha=0.01);
    ref.ami_partition_a01.(amiNew{i}) = struct('tau', tau, 'curve', curve(:,2)');
    [tau, curve] = ami_knn(y, 10);
    ref.ami_knn.(amiNew{i}) = struct('tau', tau, 'curve', curve(:,2)');
    [tau, curve] = ami_knn(y, 6, k=3);
    ref.ami_knn_k3.(amiNew{i}) = struct('tau', tau, 'curve', curve(:,2)');
    [tau, curve] = ami_copula(y, 30);
    ref.ami_copula.(amiNew{i}) = struct('tau', tau, 'curve', curve(:,2)');
end
y = readmatrix(fullfile(fx, 'ar1_phi70_512.csv'));
[tau, curve] = ami_kde(y, 8, bandwidth="adaptive");
ref.ami_kde_adaptive.ar1_phi70_512 = struct('tau', tau, 'curve', curve(:,2)');
y = readmatrix(fullfile(fx, 'lorenz_2048.csv'));
[tau, curve] = ami_kde(y, 10, bandwidth="adaptive");
ref.ami_kde_adaptive.lorenz_2048 = struct('tau', tau, 'curve', curve(:,2)');

% ---- Entropy family.
a = readmatrix(fullfile(fx, 'ar1_phi70_512.csv'));
h = readmatrix(fullfile(fx, 'henon_512.csv'));
ref.ent_xsamp.ar1_henon = [ent_xsamp(a, h, 2, 0.2, 0), ent_xsamp(a, h, 2, 0.2, 1), ...
                           ent_xsamp(a, h, 2, 0.2, 2), ent_xsamp(a, h, 1, 0.3, 2)];
ref.ent_xap.ar1_henon = [ent_xap(a, h, 2, 0.2, 1), ent_xap(a, h, 1, 0.3, 1), ...
                         ent_xap(a, h, 2, 0.25, 2)];
b = double(a > median(a));
ref.ent_symbolic.ar1_median = [ent_symbolic(b, 2), ent_symbolic(b, 3), ent_symbolic(b, 4)];
y = readmatrix(fullfile(fx, 'lorenz_2048.csv'));
Y = psr(y(1:400), 10, 3);
D = sqrt((Y(:,1) - Y(:,1)').^2 + (Y(:,2) - Y(:,2)').^2 + (Y(:,3) - Y(:,3)').^2);
ref.ent_weighted.lorenz_400 = ent_weighted(exp(-D / std(y)));
% ---- Embedding: FNN (k-d tree search, exercised at N > 512) and
% correlation dimension.
fnnCases = {'lorenz_2048', 10, 8; 'henon_512', 1, 5; 'white_4096', 1, 6; ...
            'pseudoperiodic_300', 3, 6};
for i = 1:size(fnnCases, 1)
    y = readmatrix(fullfile(fx, [fnnCases{i,1} '.csv']));
    [d0, e0] = fnn(y, fnnCases{i,2}, fnnCases{i,3}, 15, 2, 0);
    [d1, e1] = fnn(y, fnnCases{i,2}, fnnCases{i,3}, 10, 2, 1);
    ref.fnn.(fnnCases{i,1}) = struct('delay', fnnCases{i,2}, 'maxdim', fnnCases{i,3}, ...
        'dim_speed0', d0, 'dE_speed0', e0', 'dim_speed1', d1, 'dE_speed1', e1');
end
cdCases = {'lorenz_2048', 10, 3; 'henon_512', 1, 2; 'ar1_phi70_512', 1, 3};
for i = 1:size(cdCases, 1)
    y = readmatrix(fullfile(fx, [cdCases{i,1} '.csv']));
    [cod, ex] = corr_dim(y, cdCases{i,2}, cdCases{i,3}, 0);
    ref.corr_dim.(cdCases{i,1}) = struct('delay', cdCases{i,2}, 'dim', cdCases{i,3}, ...
        'CoD', cod, 'logEps', ex.logEps', 'logC', ex.logC', 'idx', ex.idx', ...
        'r2', ex.r2, 'bounds', ex.bounds);
end

% ---- RQA family. RP is summarised by its row and column sums.
rqaPack = @(RP, R) struct('rowsum', sum(RP, 2)', 'colsum', sum(RP, 1), ...
    'REC', R.REC, 'RADIUS', R.RADIUS, 'Size', R.Size, 'DET', R.DET, 'MeanL', R.MeanL, ...
    'MaxL', R.MaxL, 'EntrL', R.EntrL, 'LAM', R.LAM, 'MeanV', R.MeanV, 'MaxV', R.MaxV, ...
    'EntrV', R.EntrV, 'EntrW', localField(R, 'EntrW'));
lz = readmatrix(fullfile(fx, 'lorenz_2048.csv'));
hn = readmatrix(fullfile(fx, 'henon_512.csv'));
ar = readmatrix(fullfile(fx, 'ar1_phi70_512.csv'));
[RP, R] = rqa(lz(1:600), 10, 3, "rec", 2.5);
ref.rqa.lorenz_rec = rqaPack(RP, R);
[RP, R] = rqa(lz(1:600), 10, 3, "rad", 0.3, norm="euc", theiler=3);
ref.rqa.lorenz_rad_euc_theiler = rqaPack(RP, R);
[RP, R] = rqa(hn, 1, 2, "rec", 5, norm="max");
ref.rqa.henon_rec_max = rqaPack(RP, R);
[RP, R] = rqa(ar, 1, 1, "rad", 0.5, norm="min", dmin=3, vmin=3, zscore=0);
ref.rqa.ar1_rad_min = rqaPack(RP, R);
[RP, R] = crqa([ar(1:400), hn(1:400)], 1, 2, "rec", 3);
ref.crqa.ar1_henon_rec = rqaPack(RP, R);
[RP, R] = crqa([lz(1:400), lz(6:405)], 4, 3, "rad", 0.4, norm="euc");
ref.crqa.lorenz_rad_euc = rqaPack(RP, R);
[RP, R] = jrqa([lz(1:300), lz(11:310)], [5 5], [2 2], "rec", 5);
ref.jrqa.lorenz_rec = rqaPack(RP, R);
[RP, R] = jrqa([ar(1:300), hn(1:300)], [1 1], [1 1], "rad", 0.5, norm="max");
ref.jrqa.ar1_henon_rad = rqaPack(RP, R);
[RP, R] = mdrqa([lz(1:400), hn(1:400)], 1, 1, "rec", 4);
ref.mdrqa.lorenz_henon_rec = rqaPack(RP, R);
[RP, R] = mdrqa([lz(1:400), hn(1:400)], 2, 2, "rad", 0.5, norm="max");
ref.mdrqa.lorenz_henon_rad = rqaPack(RP, R);

% ---- Lyapunov: both estimators and the wrapper.
lz = readmatrix(fullfile(fx, 'lorenz_2048.csv'));
hn = readmatrix(fullfile(fx, 'henon_512.csv'));
o = lye_r(lz, 100, 10, 3);
ref.lye_r.lorenz = struct('pairs', o(:,2)', 'divergence', o(:,3)');
o = lye_r(lz, 100, 10, 3, 'theiler', 25);
ref.lye_r.lorenz_theiler25 = struct('pairs', o(:,2)', 'divergence', o(:,3)');
[lyS, lyL, o] = lye_r(lz, 100, 10, 3, [0 0.4 1 3], 0.75, 0);
ref.lye_r.lorenz_slopes = struct('LyES', lyS, 'LyEL', lyL);
o = lye_r(hn, 1, 1, 2);
ref.lye_r.henon = struct('pairs', o(:,2)', 'divergence', o(:,3)');
[o, L] = lye_w(lz, 100, 10, 3, 10);
ref.lye_w.lorenz_e10 = struct('out', reshape(o', 1, []), 'ncol', size(o,2), 'LyE', L);
[o, L] = lye_w(hn, 1, 1, 2, 1);
ref.lye_w.henon_e1 = struct('out', reshape(o', 1, []), 'ncol', size(o,2), 'LyE', L);
[o, L] = lye_w(lz(1:1500), 100, 8, 4, 5, 2, 20*pi/180, 1);
ref.lye_w.lorenz_custom = struct('out', reshape(o', 1, []), 'ncol', size(o,2), 'LyE', L);
[lam, ex] = lyapunov(lz, 100, delay=10, dim=3);
ref.lyapunov.lorenz_rosenstein = struct('lambda', lam, 'scalingRegion', ex.scalingRegion', ...
    'fitR2', ex.fitR2);
[lam, ex] = lyapunov(hn, 1, delay=1, dim=2);
ref.lyapunov.henon_rosenstein = struct('lambda', lam, 'scalingRegion', ex.scalingRegion', ...
    'fitR2', ex.fitR2);
lam = lyapunov(lz, 100, algorithm="wolf", delay=10, dim=3, evolve=10);
ref.lyapunov.lorenz_wolf = struct('lambda', lam);

% ---- Continuous relative phase (Signal Processing Toolbox for hilbert).
if exist('hilbert', 'file')
    lz = readmatrix(fullfile(fx, 'lorenz_2048.csv'));
    [crp, crpH] = rel_phase_cont(lz(1:800), lz(21:820), 100);
    ref.rel_phase_cont.lorenz = struct('crp', crp', 'crpH', crpH');
end

% ---- chaos_library: maps are exact arithmetic; flows go through ode45 or
% ode23, compared over a short horizon before chaos amplifies rounding.
[~, y] = chaos_library('Hennon', 500, [0.1 0.1], [1.4 0.3]);
ref.chaos_library.henon = struct('y', reshape(y', 1, []));
[~, y] = chaos_library('Logistic', 300, 0.2, 3.9);
ref.chaos_library.logistic = struct('y', y');
flows = {'Lorenz', [1 1 1], [10 28 8/3]; 'Rossler', [1 1 0], [0.2 0.2 5.7]; ...
         'Aizawa', [0.1 0 0], [0.7 0.95 0.6 3.5 0.25 0.1]; ...
         'NoseHoover', [0 5 0], 1; 'Arneodo', [0.2 0.2 -0.75], [-5.5 3.5 -1]};
for i = 1:size(flows, 1)
    [~, y] = chaos_library(flows{i,1}, 0:0.01:2, flows{i,2}, flows{i,3});
    ref.chaos_library.(lower(flows{i,1})) = struct('y', reshape(y', 1, []));
end

msCases = {'ar1_phi70_512', 6, 2, 0.15; 'henon_512', 4, 2, 0.2; 'pseudoperiodic_300', 3, 1, 0.2};
for i = 1:size(msCases, 1)
    y = readmatrix(fullfile(fx, [msCases{i,1} '.csv']));
    [rc, c, m, f, g] = ent_ms_plus(y, msCases{i,2}, msCases{i,3}, msCases{i,4});
    ref.ent_ms_plus.(msCases{i,1}) = struct('scale', msCases{i,2}, 'dim', msCases{i,3}, ...
        'radius', msCases{i,4}, 'rcmse', rc', 'cmse', c', 'mse', m', 'msfe', f', 'gmse', g');
end

ref.generated_by = version;
out = fullfile(fx, 'matlab_reference.json');
fid = fopen(out, 'w');
fprintf(fid, '%s', jsonencode(ref, 'PrettyPrint', true));
fclose(fid);
fprintf('wrote %s\n', out);
end

function v = localFirstMinAmi(y)
try
    a = ami_thomas(y, 20);
    v = a(1, 1);
catch ME
    v = NaN;
    fprintf('AMI_Thomas unavailable for reference: %s\n', ME.message);
end
end

function v = localField(S, name)
if isfield(S, name)
    v = S.(name);
else
    v = NaN;
end
end
