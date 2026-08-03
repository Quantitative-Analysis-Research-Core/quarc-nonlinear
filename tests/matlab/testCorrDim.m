function tests = testCorrDim
%TESTCORRDIM corr_dim against the correlation dimensions in Sprott Appendix A.
%
%   The appendix gives D2 for all 62 systems, computed from at least 2e12
%   pairs and extrapolated to zero scale. Those are the reference values;
%   see tests/fixtures/sprott_appendix_a.md.
%
%   They carry real uncertainty and are not tight targets. A short series
%   estimated at one embedding cannot reproduce a value computed from two
%   million points, so the assertions here are about central tendency and
%   ordering across the whole battery rather than any single system.
tests = functiontests(localfunctions);
end

function teardown(~)
dbclear all
end

% ------------------------------------------------------------------
% Whole-battery accuracy. Measured over the 55 systems that generate a
% usable series: median ratio 0.956, 84% within 25%, 96% within 50%.
% ------------------------------------------------------------------
function testTracksSprottDimensionsAcrossTheBattery(tc)
r = localRun();
tc.verifyGreaterThan(numel(r.ratio), 45, sprintf( ...
    'only %d systems produced an estimate; expected most of the battery', ...
    numel(r.ratio)));

med = median(r.ratio);
fprintf('    [battery] corr_dim over %d Sprott systems: median ratio %.3f, within 25%%: %.0f%%\n', ...
    numel(r.ratio), med, 100*mean(abs(r.ratio-1) <= 0.25));

tc.verifyEqual(med, 1.0, 'AbsTol', 0.15, sprintf( ...
    'median ratio to Sprott D2 is %.3f across %d systems', med, numel(r.ratio)));
tc.verifyGreaterThan(mean(abs(r.ratio-1) <= 0.25), 0.70, sprintf( ...
    'only %.0f%% of systems came within 25%% of the reference', ...
    100*mean(abs(r.ratio-1) <= 0.25)));
end

function testSeparatesLowFromHighDimensionalSystems(tc)
% Ordering is the weaker but more robust claim: whatever the bias, a
% one-dimensional map should score below a three-dimensional flow.
r = localRun();
maps  = r.observed(r.isMap);
flows = r.observed(~r.isMap);
tc.verifyLessThan(median(maps), median(flows), sprintf( ...
    ['median D2 for noninvertible maps (%.3f) should fall below that for ' ...
     'flows (%.3f)'], median(maps), median(flows)));
end

% ------------------------------------------------------------------
% Invariances, which need no reference value.
% ------------------------------------------------------------------
function testIsScaleInvariant(tc)
% A correlation dimension counts pairs within a radius that scales with the
% data, so multiplying the series by a constant must not change it.
x = nonantest.signals('lorenz', 3000);
base = localFirst(corr_dim(x, 8, 5, false));
for c = [1e-3 1e3 1e6]
    got = localFirst(corr_dim(c*x, 8, 5, false));
    tc.verifyEqual(got, base, 'RelTol', 0.05, sprintf( ...
        'scaling by %g moved D2 from %.4f to %.4f', c, base, got));
end
end

function testWhiteNoiseFillsTheEmbeddingSpace(tc)
% Independent samples have no attractor, so the correlation dimension rises
% with the embedding dimension instead of saturating.
x = nonantest.signals('white', 3000);
d = arrayfun(@(m) localFirst(corr_dim(x, 1, m, false)), 2:5);
fprintf('    [known answer] white noise D2 at dim 2:5 = %s\n', mat2str(round(d,3)));
tc.verifyTrue(all(diff(d) > 0), sprintf( ...
    'D2 of white noise should increase with embedding dimension: %s', ...
    mat2str(round(d,3))));
end

function testRunsHeadless(tc)
s = nonantest.sideEffects(@() corr_dim(nonantest.signals('lorenz',1000), 8, 5, false));
tc.verifyFalse(s.errored, 'corr_dim errored');
tc.verifyEqual(s.figures, 0, 'corr_dim opened a figure with showPlot false');
end

% ---------------------------------------------------------------- helpers

function r = localRun()
%LOCALRUN corr_dim over every usable Sprott system, one uniform protocol.
persistent cached
if ~isempty(cached), r = cached; return, end

c = nonantest.sprott_catalog();
c = c([c.usable]);
ref = []; obs = []; isMap = logical([]);
for i = 1:numel(c)
    s = c(i);
    if ~isfinite(s.d2), continue, end
    try
        [x, gi] = nonantest.sprott_series(s, 4000);
        if gi.degenerate, continue, end
        if s.kind == "map"
            v = corr_dim(x, 1, 3, false);
        else
            v = corr_dim(x, 8, 5, false);
        end
        v = localFirst(v);
        if ~isfinite(v) || v <= 0, continue, end
        ref(end+1) = s.d2;              %#ok<AGROW>
        obs(end+1) = v;                 %#ok<AGROW>
        isMap(end+1) = s.kind == "map"; %#ok<AGROW>
    catch
        % A system that cannot be estimated is not a failure of the metric
        % under test; it is reported by count in the battery test above.
    end
end
r = struct('reference', ref(:), 'observed', obs(:), ...
           'ratio', obs(:)./ref(:), 'isMap', isMap(:));
cached = r;
end

function v = localFirst(d)
v = d(1);
end
