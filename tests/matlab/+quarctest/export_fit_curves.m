function paths = export_fit_curves(opts)
%EXPORT_FIT_CURVES Write the curves both scaling-region fits are taken from.
%
%   paths = quarctest.export_fit_curves()
%   paths = quarctest.export_fit_curves(Systems=["lorenz","lcg"], N=4000)
%
%   Writes two long-format CSVs into tests/reports:
%
%     curves_corr_dim.csv    ln(epsilon), ln(C(epsilon)), and whether each
%                            point fell inside the fitted band
%     curves_lyapunov.csv    the Rosenstein divergence curve, sample index
%                            against mean log divergence, and whether each
%                            point fell inside the fitted window
%
%   WHY THESE EXIST. The characterization report gives a slope per system and,
%   now, an R^2 for the region it was fitted over. Neither shows the SHAPE of
%   the curve, and a correlation dimension or a Lyapunov exponent is only
%   meaningful if the stretch it was fitted to is actually straight. These
%   files carry the curves themselves so the claim can be looked at rather
%   than summarised.
%
%   ONE REALIZATION PER SYSTEM, AND IT IS SPROTT'S. The published initial
%   condition is realization 1 of every ensemble, so the curve written here is
%   a member of the distribution the report summarises rather than a separate
%   run. The protocol matches quarctest.characterize exactly: same decimation,
%   same embedding policy, same series length.
%
%   Options
%     Systems  names to export. Default spans the range of fit quality seen
%              across the catalogue rather than only the well-behaved cases.
%     N        samples per series. Default 4000, matching the full run.
%
%   See also QUARCTEST.CHARACTERIZE, QUARCTEST.METRIC_POLICY.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    opts.Systems (1,:) string = ["logistic", "henon", "lorenz", "rossler", ...
                                 "lcg", "duffing_two_well", "henon_area", "thomas"]
    opts.N       (1,1) double {mustBePositive, mustBeInteger} = 4000
end

quarctest.require_library();
c = quarctest.sprott_catalog();
c = c([c.usable]);

cdRows = {};
lyRows = {};

for k = 1:numel(opts.Systems)
    sys = c([c.name] == opts.Systems(k));
    if isempty(sys)
        error('quarctest:exportFitCurves:unknownSystem', ...
              '%s is not in the usable catalogue', opts.Systems(k));
    end
    E = quarctest.embed_policy(sys);

    [~, gi0] = quarctest.sprott_series(sys, min(opts.N, 2048));
    decim = gi0.decim;
    if sys.kind == "map"
        [x, gi] = quarctest.sprott_series(sys, opts.N);
    else
        [x, gi] = quarctest.sprott_series(sys, opts.N, Decim=decim);
    end
    if gi.degenerate, continue, end

    if E.delayRule == "fixed", delay = E.delayValue; else, delay = ami(x, 100); end
    if ~isfinite(delay) || delay < 1, delay = 1; end
    if E.dimRule == "fixed"
        dim = E.dimValue;
    else
        dim = fnn(x, delay, 10, 15, 2, 1);
        if ~isfinite(dim) || dim < E.stateDim, dim = E.stateDim + 1; end
    end

    % ---- correlation sum curve
    try
        [~, cex] = corr_dim(x, delay, dim, false);
        inFit = false(numel(cex.logEps), 1);
        inFit(cex.idx) = true;
        for j = 1:numel(cex.logEps)
            cdRows{end+1} = {sys.name, sys.category, cex.logEps(j), cex.logC(j), ...
                             double(inFit(j)), cex.r2, cex.slope, sys.d2}; %#ok<AGROW>
        end
    catch
    end

    % ---- Rosenstein divergence curve
    try
        [~, lex] = lyapunov(x, gi.fs, algorithm="rosenstein", delay=delay, dim=dim);
        d = lex.divergence(:);
        inFit = false(numel(d), 1);
        inFit(lex.scalingRegion) = true;
        for j = 1:numel(d)
            lyRows{end+1} = {sys.name, sys.category, j, d(j), ...
                             double(inFit(j)), lex.fitR2, gi.fs, sys.lambda}; %#ok<AGROW>
        end
    catch
    end
end

here = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % tests/
outDir = fullfile(here, 'reports');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cdPath = fullfile(outDir, 'curves_corr_dim.csv');
lyPath = fullfile(outDir, 'curves_lyapunov.csv');

T1 = cell2table(vertcat(cdRows{:}), 'VariableNames', ...
    {'system','category','logEps','logC','inFit','fitR2','slope','referenceD2'});
T2 = cell2table(vertcat(lyRows{:}), 'VariableNames', ...
    {'system','category','sample','divergence','inFit','fitR2','fs','referenceLambda'});

writetable(T1, cdPath);
writetable(T2, lyPath);

paths = struct('corrDim', string(cdPath), 'lyapunov', string(lyPath));
fprintf('wrote %s (%d rows)\n', cdPath, height(T1));
fprintf('wrote %s (%d rows)\n', lyPath, height(T2));
end
