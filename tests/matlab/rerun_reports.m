function rerun_reports(stage)
%RERUN_REPORTS Regenerate the committed characterization reports, restartably.
%
%   matlab -batch "addpath('tests/matlab'); rerun_reports('dysts')"
%
%   Stages, in the order tests/rerun_reports.sh runs them:
%
%     fast        characterize(R=10, N=2000)          characterization_fast.*
%     curves      export_fit_curves()                 curves_*.csv
%     evolve      evolve_sweep()                      evolve_sweep.csv
%     full        characterize(R=1000, N=4000)        characterization_full.*
%     dysts       characterize_series(R=100)          characterization_dysts.*
%     noise       noise_sweep(Family="lyapunov")      noise_sweep.csv
%     noiserqa    noise_sweep(Family="rqa")           noise_sweep_rqa/chunk_NN.csv
%
%   The arguments are the ones recorded in each report's header. The long
%   stages run in chunks of a few systems, one .mat per chunk under
%   tests/artifacts/rerun/<stage>/, and a chunk already on disk is skipped,
%   so a stage killed partway resumes where it stopped. Chunking does not
%   change any value: every random draw is seeded from the system's name
%   (quarctest.ensemble_ics), not from its position in the run.
%
%   noise_sweep is the exception: it seeds each system's noise from the
%   system's position within the call, so its stages reproduce the original
%   call pattern exactly. The Lyapunov pass is one call, as committed; the
%   RQA pass keeps the original chunks (manifest systems in alphabetical
%   order, six per chunk), merged by python/merge_noise_sweep_chunks.py, as
%   in the 28-30 September run.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(here);
addpath(genpath(fullfile(repo, 'matlab')));
reports = fullfile(repo, 'tests', 'reports');
work = fullfile(repo, 'tests', 'artifacts', 'rerun', stage);
if ~exist(work, 'dir'), mkdir(work); end
dbclear all

switch stage
    case 'fast'
        [S, p] = quarctest.characterize(R=10, N=2000);
        quarctest.write_characterization(S, p, Tag='fast');

    case 'curves'
        quarctest.export_fit_curves();

    case 'evolve'
        T = quarctest.evolve_sweep();
        writetable(T, fullfile(reports, 'evolve_sweep.csv'));

    case 'full'
        C = quarctest.sprott_catalog();
        C = C([C.usable]);                 % as characterize filters them
        [S, p] = localChunks(work, string({C.name}), 5, ...
            @(names) localPair(@() quarctest.characterize(R=1000, N=4000, Systems=names)));
        quarctest.write_characterization(S, p, Tag='full');

    case 'dysts'
        [S, p] = localChunks(work, localManifestSystems(reports), 6, ...
            @(names) localPair(@() quarctest.characterize_series(R=100, Systems=names)));
        quarctest.write_characterization(S, p, Tag='dysts');

    case 'noise'
        % One call, not chunks: noise_sweep seeds each system's noise from
        % its position in the call (Seed + i*100000), so chunking would draw
        % different noise from the committed single-call run. ~7 h.
        [S, p] = quarctest.noise_sweep(Family="lyapunov", ...
            Checkpoint=fullfile(work, 'noise_sweep.ckpt.mat'));
        writetable(S, fullfile(reports, 'noise_sweep.csv'));
        writetable(p, fullfile(reports, 'noise_sweep_per_realization.csv'));

    case 'noiserqa'
        names = sort(localManifestSystems(reports));
        outDir = fullfile(reports, 'noise_sweep_rqa');
        perDir = fullfile(repo, 'tests', 'artifacts', 'noise_sweep_rqa');
        if ~exist(perDir, 'dir'), mkdir(perDir); end
        nChunk = ceil(numel(names) / 6);
        for k = 1:nChunk
            done = fullfile(work, sprintf('chunk_%02d.done', k));
            if isfile(done), continue, end
            sel = names((k-1)*6 + 1 : min(k*6, numel(names)));
            fprintf('[noiserqa] chunk %d/%d: %s\n', k, nChunk, strjoin(sel, ', '));
            [S, p] = quarctest.noise_sweep(Family="rqa", Systems=sel);
            writetable(S, fullfile(outDir, sprintf('chunk_%02d.csv', k)));
            writetable(p, fullfile(perDir, sprintf('chunk_%02d.csv', k)));
            fclose(fopen(done, 'w'));
        end

    otherwise
        error('rerun_reports:unknownStage', 'Unknown stage "%s".', stage);
end
fprintf('[%s] done\n', stage);
end

% ---------------------------------------------------------------- helpers

function [S, p] = localChunks(work, names, per, run)
%LOCALCHUNKS Run RUN over NAMES in chunks, caching each chunk's tables.
nChunk = ceil(numel(names) / per);
Sall = cell(nChunk, 1);
pall = cell(nChunk, 1);
for k = 1:nChunk
    f = fullfile(work, sprintf('chunk_%02d.mat', k));
    if isfile(f)
        c = load(f);
    else
        sel = names((k-1)*per + 1 : min(k*per, numel(names)));
        fprintf('[chunk %d/%d] %s\n', k, nChunk, strjoin(sel, ', '));
        [c.S, c.p] = run(sel);
        save(f, '-struct', 'c');
    end
    Sall{k} = c.S;
    pall{k} = c.p;
end
S = vertcat(Sall{:});
p = vertcat(pall{:});
end

function [S, p] = localPair(fn)
[S, p] = fn();
end

function names = localManifestSystems(reports)
m = jsondecode(fileread(fullfile(reports, 'dysts', 'manifest.json')));
names = string({m.systems.system});
end
