function tests = testEntCross
%TESTENTCROSS Contract tests for ent_ap, ent_xap, ent_xsamp and ent_weighted.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
tc.TestData.white = nonantest.signals('white', 1000);
tc.TestData.ar1   = nonantest.signals('ar1', 1000, 0.8);
tc.TestData.sine  = nonantest.signals('sine', 1000, 25);
end

function teardown(~)
dbclear all
end

% ------------------------------------------------------------------
% ent_xsamp's third normalisation option assigned nothing.
%
% The branch printed "These data will not be normalized" and returned without
% setting xn, yn or r, so the comparison loop below it failed on an undefined
% variable. The option is documented and was unusable.
% ------------------------------------------------------------------
function testAllNormalisationOptionsWork(tc)
x = tc.TestData.white; y = tc.TestData.ar1;
got = nan(1,3);
for n = 1:3
    s = nonantest.sideEffects(@() ent_xsamp(x, y, 2, 0.2, n));
    tc.verifyFalse(s.errored, sprintf( ...
        'ent_xsamp with norm=%d errored: %s', n, localMsg(s)));
    if ~s.errored, got(n) = s.value; end
end
tc.verifyTrue(all(isfinite(got)), sprintf( ...
    'not every normalisation produced a finite value: %s', mat2str(round(got,4))));
end

% ------------------------------------------------------------------
% Ordering by regularity. These hold for any sensible entropy.
% ------------------------------------------------------------------
function testApproximateEntropyOrdersByRegularity(tc)
white = ent_ap(tc.TestData.white, 2, 0.2);
ar1   = ent_ap(tc.TestData.ar1,   2, 0.2);
sine  = ent_ap(tc.TestData.sine,  2, 0.2);
fprintf('    [known answer] ent_ap: white %.4f, ar1 %.4f, sine %.4f\n', white, ar1, sine);
tc.verifyGreaterThan(white, ar1, 'white noise should exceed AR(1)');
tc.verifyGreaterThan(ar1, sine, 'AR(1) should exceed a pure sine');
tc.verifyLessThan(sine, 0.5, sprintf( ...
    'a pure sine is highly regular; ent_ap gave %.4f', sine));
end

function testCrossEntropiesRiseWithDissimilarity(tc)
% Two independent series share fewer template matches than a series and a
% lightly perturbed copy of itself.
x = tc.TestData.ar1;
near = x + 0.05*std(x)*nonantest.signals('white', numel(x));
far  = tc.TestData.white;
for f = {@ent_xsamp, @ent_xap}
    fn = f{1};
    tc.verifyLessThan(fn(x, near, 2, 0.2, 1), fn(x, far, 2, 0.2, 1), sprintf( ...
        '%s did not score a near copy below an independent series', func2str(fn)));
end
end

% ------------------------------------------------------------------
% Invariances.
% ------------------------------------------------------------------
function testCrossEntropiesAreSymmetric(tc)
% Cross entropy measures shared structure, which is a property of the pair.
x = tc.TestData.white; y = tc.TestData.ar1;
tc.verifyEqual(ent_xsamp(x,y,2,0.2,1), ent_xsamp(y,x,2,0.2,1), 'RelTol', 1e-12, ...
    'ent_xsamp must be symmetric in its two series');
tc.verifyEqual(ent_xap(x,y,2,0.2,1), ent_xap(y,x,2,0.2,1), 'RelTol', 1e-12, ...
    'ent_xap must be symmetric in its two series');
end

function testScaleInvariance(tc)
% The radius is a proportion of the spread, so a change of units must not
% change the result.
x = tc.TestData.ar1; y = tc.TestData.white;
tc.verifyEqual(ent_ap(1000*x,2,0.2), ent_ap(x,2,0.2), 'RelTol', 1e-10, ...
    'ent_ap is not scale invariant');
tc.verifyEqual(ent_xsamp(1000*x,1000*y,2,0.2,1), ent_xsamp(x,y,2,0.2,1), ...
    'RelTol', 1e-10, 'ent_xsamp is not scale invariant');
end

function testCrossApEnIncrementArgument(tc)
% k is the increment between the two template lengths compared: the function
% evaluates dim and dim+k, so k=1 is the standard cross approximate entropy.
% The header described it only as "something lag".
x = tc.TestData.white; y = tc.TestData.ar1;
v = arrayfun(@(k) ent_xap(x, y, 2, 0.2, k), 1:3);
tc.verifyTrue(all(isfinite(v)), 'ent_xap produced non-finite values');
tc.verifyTrue(all(diff(v) > 0), sprintf( ...
    ['comparing more widely separated template lengths should not reduce ' ...
     'the entropy: %s'], mat2str(round(v,4))));
end

% ------------------------------------------------------------------
% ent_weighted.
% ------------------------------------------------------------------
function testWeightedEntropyBounds(tc)
% Uniform column weights occupy a single bin and carry no information.
tc.verifyEqual(ent_weighted(ones(60)), 0, 'AbsTol', 1e-12, ...
    'a uniform weighted recurrence plot has zero entropy');
rng(3);
tc.verifyGreaterThan(ent_weighted(-abs(randn(60))), 0, ...
    'varied column weights should give positive entropy');
end

function testWeightedEntropyIsScaleInvariant(tc)
% The 49 bins span the observed range of column weights, so the measure
% depends on the shape of their distribution and not on its units.
rng(5);
W = -abs(randn(60));
base = ent_weighted(W);
for c = [10 1000]
    tc.verifyEqual(ent_weighted(c*W), base, 'RelTol', 1e-12, sprintf( ...
        'scaling the weights by %g moved the entropy from %.6f to %.6f', ...
        c, base, ent_weighted(c*W)));
end
end

function testWeightedEntropyDistinguishesConcentratedFromSpread(tc)
% Since the bins are relative to the range, what the measure sees is how the
% column weights are distributed within it: mass in one bin is low entropy,
% mass spread across bins is high.
n = 60;
concentrated = -ones(n);
concentrated(:,1) = -100;                              % one outlier column
spread = -abs(repmat(linspace(1,100,n), n, 1));        % evenly spread sums
lo = ent_weighted(concentrated);
hi = ent_weighted(spread);
fprintf('    [known answer] ent_weighted: concentrated %.4f, spread %.4f\n', lo, hi);
tc.verifyLessThan(lo, hi, sprintf( ...
    ['weights concentrated in one bin scored %.4f, evenly spread weights ' ...
     '%.4f; concentrated should be lower'], lo, hi));
end

% ------------------------------------------------------------------
% Housekeeping.
% ------------------------------------------------------------------
function testFamilyIsSilent(tc)
x = nonantest.signals('white', 400); y = nonantest.signals('ar1', 400, 0.5);
calls = { ...
    'ent_ap',        @() ent_ap(x, 2, 0.2); ...
    'ent_xap',       @() ent_xap(x, y, 2, 0.2, 1); ...
    'ent_xsamp n=1', @() ent_xsamp(x, y, 2, 0.2, 1); ...
    'ent_xsamp n=3', @() ent_xsamp(x, y, 2, 0.2, 3); ...
    'ent_weighted',  @() ent_weighted(-abs(randn(30)))};
for k = 1:size(calls,1)
    out = evalc('calls{k,2}();');
    tc.verifyEmpty(strtrim(out), sprintf( ...
        '%s wrote to the console: "%s"', calls{k,1}, strtrim(out)));
end
end

function testRunHeadless(tc)
x = nonantest.signals('white', 400); y = nonantest.signals('ar1', 400, 0.5);
calls = {@() ent_ap(x,2,0.2), @() ent_xap(x,y,2,0.2,1), ...
         @() ent_xsamp(x,y,2,0.2,1), @() ent_weighted(-abs(randn(30)))};
for k = 1:numel(calls)
    s = nonantest.sideEffects(calls{k});
    tc.verifyFalse(s.errored, sprintf('call %d errored: %s', k, localMsg(s)));
    tc.verifyEqual(s.figures, 0, sprintf('call %d opened a figure', k));
end
end

% ---------------------------------------------------------------- helper

function m = localMsg(s)
if isempty(s.err), m = ''; else, m = s.err.message; end
end
