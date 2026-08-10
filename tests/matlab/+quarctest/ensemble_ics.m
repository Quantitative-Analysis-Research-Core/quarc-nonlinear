function [ics, info] = ensemble_ics(sys, R, opts)
%ENSEMBLE_ICS I.i.d. initial conditions around a system's published value.
%
%   [ics, info] = quarctest.ensemble_ics(sys, R)
%   [ics, info] = quarctest.ensemble_ics(sys, R, Seed=S, Spread=P)
%
%   Returns a d-by-R matrix whose columns are initial conditions for the
%   system sys, and a struct recording exactly how they were drawn.
%
%   COLUMN 1 IS THE PUBLISHED INITIAL CONDITION, VERBATIM. Sprott's own x0 is
%   always a member of the ensemble, so the published trajectory is one of the
%   realizations rather than a separate special case, and every summary
%   statistic contains it.
%
%   Columns 2:R are x0 + Spread*scale.*randn, drawn i.i.d., where scale is
%   max(abs(x0), 1) per component so that a component published as 0 is
%   perturbed on the same absolute footing as one published as 6.
%
%   WHAT THE SPREAD IS FOR, AND WHAT IT IS NOT. For a dissipative system the
%   transient discarded by quarctest.sprott_series erases the initial
%   condition: every realization lands on the same attractor, at a different
%   phase. The spread therefore selects WHERE ON THE ATTRACTOR a realization
%   starts, and the resulting scatter in a metric is that estimator's sampling
%   variability over finite windows of one attractor. It is not a measure of
%   sensitivity to initial conditions, and it does not explore the basin.
%
%   The default spread is deliberately small for that reason. A multistable
%   system has more than one attractor, and a spread large enough to cross a
%   basin boundary would produce an ensemble whose scatter mixes two different
%   objects -- a bimodal metric distribution reported as a mean and an SD,
%   which is exactly the kind of number that means nothing. Conservative
%   systems have no attractor at all, so their realizations are genuinely
%   distinct orbits and their spread is a different quantity again; the report
%   marks the category so the two are not read as the same measurement.
%
%   REPRODUCIBILITY. The draw is a pure function of (Seed, system name, R,
%   Spread). info.seed and info.streamSeed record what was used, and calling
%   this function again with the same arguments returns the identical matrix.
%   The per-system stream seed is derived from the system name, so adding,
%   removing or reordering systems in the catalogue does not change the
%   initial conditions drawn for any other system.
%
%   NO TRAJECTORY IS RUN HERE. An initial condition that produces a
%   degenerate orbit is not filtered out at this stage, because silently
%   redrawing would make the ensemble size depend on the rejection rule and
%   hide how often it fires. Degeneracy is detected downstream by
%   quarctest.sprott_series and reported as a realization count below R.
%
%   Fields of info
%     seed        the Seed argument as supplied
%     streamSeed  the per-system seed actually used for the draw
%     spread      the Spread argument as supplied
%     scale       per-component perturbation scale, max(abs(x0),1)
%     published   the published initial condition (== ics(:,1))
%
%   See also QUARCTEST.SPROTT_CATALOG, QUARCTEST.SPROTT_SERIES,
%   QUARCTEST.CHARACTERIZE.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    sys (1,1) struct
    R   (1,1) double {mustBePositive, mustBeInteger}
    opts.Seed   (1,1) double = 20260810
    opts.Spread (1,1) double {mustBeNonnegative} = 0.01
end

x0 = sys.x0(:);
d  = numel(x0);

% Derive a per-system stream seed from the name so that the ensemble for one
% system does not depend on how many systems precede it in the catalogue.
streamSeed = mod(opts.Seed + localNameHash(sys.name), 2^31 - 1);

s = RandStream('mt19937ar', 'Seed', streamSeed);

scale = max(abs(x0), 1);

ics = zeros(d, R);
ics(:,1) = x0;
if R > 1
    ics(:,2:R) = x0 + opts.Spread .* scale .* randn(s, d, R-1);
end

info = struct('seed', opts.Seed, 'streamSeed', streamSeed, ...
              'spread', opts.Spread, 'scale', scale, 'published', x0);
end

% ---------------------------------------------------------------- internals

function h = localNameHash(name)
%LOCALNAMEHASH Small deterministic hash of a system name.
%   djb2, reduced mod 2^31-1. Any stable integer-valued function of the name
%   would do; what matters is that it does not depend on catalogue order and
%   does not change between MATLAB releases, which rules out string2hash-style
%   helpers and anything built on java.lang.String.hashCode.
c = double(char(name));
h = 5381;
for i = 1:numel(c)
    h = mod(h * 33 + c(i), 2^31 - 1);
end
end
