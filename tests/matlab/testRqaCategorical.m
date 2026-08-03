function tests = testRqaCategorical
%TESTRQACATEGORICAL RQA on binary and categorical series.
%   Recurrence analysis is routinely run on symbolic data -- binary event
%   sequences, categorical states, discretised signals -- where the distance
%   matrix takes only a handful of distinct values. That regime broke several
%   assumptions that hold for continuous data.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
x = nonantest.signals('lorenz', 200);
tc.TestData.continuous  = x;
tc.TestData.binary      = double(x > 0);
tc.TestData.categorical = double(discretize(x, 4));
end

function teardown(~)
dbclear all
end

% ------------------------------------------------------------------
% The distance matrix takes few distinct values, so %REC is a step function
% of radius and an arbitrary target falls between steps.
% ------------------------------------------------------------------
function testUnreachableTargetFailsInsteadOfHanging(tc)
% Searching for %REC = 2.5 on a binary series cannot succeed: %REC jumps from
% 0 to ~50 as the radius crosses zero. The search previously halved the radius
% without bound and never terminated -- measured still running after 90 s,
% with the radius down to 6e-38 and %REC stuck at 49.7.
s = nonantest.sideEffects(@() rqa(tc.TestData.binary, 1, 1, "rec", 2.5, Norm="none"));
tc.verifyTrue(s.errored, ...
    'An unreachable recurrence target must raise, not search forever.');
tc.verifyEqual(s.err.identifier, 'set_radius:targetUnreachable', ...
    sprintf('unexpected error: %s', s.err.message));
tc.verifyLessThan(s.seconds, 30, sprintf( ...
    'took %.1f s to give up; it should fail promptly', s.seconds));
end

function testExplicitRadiusWorksOnSymbolicData(tc)
% The documented alternative to a recurrence target. It must actually work,
% since the error message above tells the caller to use it.
for name = ["binary", "categorical"]
    y = tc.TestData.(name);
    [rp, r] = rqa(y, 1, 1, "rad", 0.5, Norm="none");
    tc.verifyNotEmpty(rp, sprintf('%s: empty recurrence plot', name));
    tc.verifyGreaterThan(r.REC, 0, sprintf('%s: zero recurrence', name));
    tc.verifyEqual(r.RADIUS, 0.5, sprintf( ...
        '%s: RADIUS should echo the radius supplied', name));
end
end

% ------------------------------------------------------------------
% RESULTS.RADIUS was read but never assigned in the param="rad" branch, so
% every direct-radius call failed with an undefined-variable error.
% ------------------------------------------------------------------
function testRadiusBranchAssignsRadius(tc)
x = tc.TestData.continuous;
y = [x nonantest.signals('rossler', 200)];
calls = { ...
    'rqa',   @() rqa(x, 3, 3, "rad", 1, Norm="none"); ...
    'crqa',  @() crqa(y, 3, 3, "rad", 1, Norm="none"); ...
    'mdrqa', @() mdrqa(y, 3, 3, "rad", 1, Norm="none"); ...
    'jrqa',  @() jrqa(y, [3 3], [3 3], "rad", 1, Norm="none")};
for k = 1:size(calls,1)
    s = nonantest.sideEffects(calls{k,2});
    tc.verifyFalse(s.errored, sprintf('%s with param="rad" errored: %s', ...
        calls{k,1}, localMsg(s)));
end
end

% ------------------------------------------------------------------
% ent_weighted divides the range of column weights into 49 bins. Equal
% weights give a zero bin size, and si_min:0:si_max is an empty range in
% MATLAB, so the histogram loop never ran and its accumulator was never
% created.
% ------------------------------------------------------------------
function testWeightedEntropyHandlesUniformWeights(tc)
s = nonantest.sideEffects(@() ent_weighted(ones(50)));
tc.verifyFalse(s.errored, sprintf( ...
    'ent_weighted failed on a uniform matrix: %s', localMsg(s)));
tc.verifyEqual(s.value, 0, 'AbsTol', 1e-12, ...
    'a single occupied bin carries no information, so entropy is 0');

% And it still measures something when the weights vary.
rng(1);
w = -abs(randn(50));
tc.verifyGreaterThan(ent_weighted(w), 0, ...
    'entropy should be positive when column weights differ');
end

% ------------------------------------------------------------------
% Continuous data must be unaffected by any of the above.
% ------------------------------------------------------------------
function testContinuousDataUnaffected(tc)
x = tc.TestData.continuous;
[~, r] = rqa(x, 5, 3, "rec", 2.5, Norm="none");
tc.verifyEqual(r.REC, 2.5, 'AbsTol', 0.05, sprintf( ...
    'target recurrence search returned %%REC = %.4f', r.REC));
tc.verifyGreaterThan(r.EntrW, 0, 'weighted entropy should be positive here');
end

% ------------------------------------------------------------------ helper

function m = localMsg(s)
if isempty(s.err), m = ''; else, m = s.err.message; end
end
