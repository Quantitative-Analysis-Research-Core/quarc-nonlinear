function T = evolve_sweep(opts)
%EVOLVE_SWEEP Wolf's exponent against its renormalisation interval.
%
%   T = quarctest.evolve_sweep()
%   T = quarctest.evolve_sweep(R=200, Evolve=[1 2 3 5 8 10 15 20 30])
%
%   Returns a long table: one row per system and evolve value, with the
%   ensemble median exponent, its spread, and the ratio to Sprott's published
%   value.
%
%   THE QUESTION. lye_w advances its reference point by EVOLVE samples and
%   replaces the neighbour at every step, so EVOLVE is the time the pair is
%   allowed to separate before renormalisation. Wolf's method assumes the pair
%   stays in the linear regime over that interval. Whether it does depends on
%   the system's own divergence rate, and lyapunov.m applies one default of 10
%   samples to every system in the catalogue.
%
%   Measured, that default is not a constant in any meaningful unit. On flows
%   the protocol decimates to about 40 samples per dominant period, so 10
%   samples is roughly a quarter turn of the attractor. On maps there is no
%   decimation and 10 samples is 10 iterations: on the logistic map, with
%   lambda = ln 2 per iteration, a pair separates by a factor of 2^10 = 1024
%   before it is replaced, which is far outside any linear regime. If that is
%   what drives Wolf's erratic behaviour on maps -- 0.33x on Arnold's cat,
%   0.59x on Ikeda -- then shrinking EVOLVE should move those ratios towards 1
%   and leave the flows comparatively unchanged.
%
%   THE SERIES IS GENERATED ONCE PER REALIZATION AND REUSED ACROSS EVOLVE
%   VALUES, so the comparison is within-realization: differences between evolve
%   values cannot come from having drawn different data. The ensemble uses the
%   same initial conditions as quarctest.characterize at the same seed.
%
%   Options
%     Systems  names to sweep. Default: every usable map, plus flows for
%              contrast, because a result that only looked at maps could not
%              distinguish "maps are different" from "the default is wrong".
%     Evolve   renormalisation intervals to try, in samples.
%     R        realizations per system. Default 200.
%     N        samples per realization. Default 4000.
%
%   See also QUARCTEST.CHARACTERIZE, LYAPUNOV.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    opts.Systems (1,:) string = "maps+flows"
    opts.Evolve  (1,:) double = [1 2 3 5 8 10 15 20 30]
    opts.R       (1,1) double {mustBePositive, mustBeInteger} = 200
    opts.N       (1,1) double {mustBePositive, mustBeInteger} = 4000
    opts.Seed    (1,1) double = 20260810
    opts.Verbose (1,1) logical = true
end

quarctest.require_library();
c = quarctest.sprott_catalog();
c = c([c.usable]);

if isscalar(opts.Systems) && opts.Systems == "maps+flows"
    isMap = [c.kind] == "map";
    flows = ["lorenz","rossler","chen","thomas","ueda","driven_vdp", ...
             "halvorsen","burke_shaw"];
    keep = isMap | ismember([c.name], flows);
    c = c(keep);
elseif ~(isscalar(opts.Systems) && opts.Systems == "all")
    c = c(ismember([c.name], opts.Systems));
end

ev = opts.Evolve;
rows = {};
t0 = tic;

for i = 1:numel(c)
    sys = c(i);
    E = quarctest.embed_policy(sys);
    [~, gi0] = quarctest.sprott_series(sys, min(opts.N, 2048));
    decim = gi0.decim;
    ics = quarctest.ensemble_ics(sys, opts.R, Seed=opts.Seed);

    L = nan(opts.R, numel(ev));
    N = opts.N;
    parfor r = 1:opts.R
        L(r,:) = localOne(sys, ics(:,r), N, decim, E, ev); %#ok<PFBNS>
    end

    for k = 1:numel(ev)
        v = L(:,k); v = v(isfinite(v));
        if isempty(v)
            med = NaN; mad_ = NaN; nn = 0;
        else
            med = median(v); mad_ = median(abs(v - med)); nn = numel(v);
        end
        rows{end+1} = {sys.name, string(sys.kind), sys.category, ev(k), nn, ...
                       med, mad_, sys.lambda, med/sys.lambda}; %#ok<AGROW>
    end
    if opts.Verbose
        best = NaN; bestEv = NaN;
        for k = 1:numel(ev)
            rt = rows{end-numel(ev)+k}{9};
            if isfinite(rt) && (~isfinite(best) || abs(log(rt)) < abs(log(best)))
                best = rt; bestEv = ev(k);
            end
        end
        d10 = rows{end-numel(ev)+find(ev==10,1)}{9};
        fprintf('%-22s %-5s ratio@10 %6.3f | best %6.3f at evolve=%2d | %.0fs\n', ...
                sys.name, sys.kind, d10, best, bestEv, toc(t0));
    end
end

T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'system','kind','category','evolve','n','median','mad','reference','ratio'});
end

% ---------------------------------------------------------------- internals

function out = localOne(sys, ic, N, decim, E, ev)
out = nan(1, numel(ev));
s = sys; s.x0 = ic(:);
try
    if s.kind == "map"
        [x, gi] = quarctest.sprott_series(s, N);
    else
        [x, gi] = quarctest.sprott_series(s, N, Decim=decim);
    end
catch
    return
end
if gi.degenerate, return, end

if E.delayRule == "fixed"
    delay = E.delayValue;
else
    try, delay = ami(x, 100); catch, delay = 1; end
    if ~isfinite(delay) || delay < 1, delay = 1; end
end
if E.dimRule == "fixed"
    dim = E.dimValue;
else
    try, dim = fnn(x, delay, 10, 15, 2, 1); catch, dim = NaN; end
    if ~isfinite(dim) || dim < E.stateDim, dim = E.stateDim + 1; end
end

for k = 1:numel(ev)
    try
        out(k) = lyapunov(x, gi.fs, algorithm="wolf", delay=delay, dim=dim, ...
                          evolve=ev(k));
    catch
    end
end
end
