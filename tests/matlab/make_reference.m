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
