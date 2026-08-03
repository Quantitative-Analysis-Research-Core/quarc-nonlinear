function tests = testEnt
%TESTENT Contract tests for the ent wrapper and the entropy family.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
tc.TestData.white = nonantest.signals('white', 2000);
tc.TestData.ar1   = nonantest.signals('ar1', 2000, 0.7);
tc.TestData.sine  = nonantest.signals('sine', 2000, 25);
end

function teardown(~)
dbclear all
end

% ------------------------------------------------------------------
% Permutation entropy compared a 1-by-dim row against a dim-by-1 column, so
% the difference broadcast to a dim-by-dim matrix and the `if` required every
% element to be zero. No ordinal pattern was ever counted, the histogram
% stayed empty, and the function returned exactly 0 for every input at every
% dimension. It worked only when handed a row vector, which no other function
% in the library expects.
% ------------------------------------------------------------------
function testPermutationEntropyIsNotIdenticallyZero(tc)
for dim = 2:5
    v = ent(tc.TestData.white, algorithm="permutation", dim=dim);
    tc.verifyGreaterThan(v, 0.5, sprintf( ...
        ['permutation entropy of white noise at dim=%d is %.4f. It should ' ...
         'approach log2(dim!) = %.4f, not 0.'], dim, v, log2(factorial(dim))));
end
end

function testPermutationEntropyIsOrientationIndependent(tc)
x = tc.TestData.ar1;
tc.verifyEqual(ent_permu(x(:), 3, 1), ent_permu(x(:).', 3, 1), 'RelTol', 1e-12, ...
    'a row and a column holding the same series must give the same entropy');
end

function testPermutationEntropyApproachesItsMaximum(tc)
% Independent samples visit all dim! orderings with equal probability, so the
% entropy approaches log2(dim!). Measured 2.5839 against 2.5850 at dim=3.
for dim = 3:5
    v = ent(tc.TestData.white, algorithm="permutation", dim=dim);
    mx = log2(factorial(dim));
    tc.verifyLessThanOrEqual(v, mx + 1e-9, ...
        sprintf('dim=%d: %.4f exceeds the maximum %.4f', dim, v, mx));
    tc.verifyGreaterThan(v, 0.90*mx, ...
        sprintf('dim=%d: %.4f is far below the maximum %.4f', dim, v, mx));
end
end

function testPermutationEntropyOrdersSignalsByRegularity(tc)
dim = 4;
white = ent(tc.TestData.white, algorithm="permutation", dim=dim);
ar1   = ent(tc.TestData.ar1,   algorithm="permutation", dim=dim);
sine  = ent(tc.TestData.sine,  algorithm="permutation", dim=dim);
fprintf('    [known answer] permutation entropy dim=4: white %.4f, ar1 %.4f, sine %.4f (max %.4f)\n', ...
    white, ar1, sine, log2(factorial(dim)));
tc.verifyGreaterThan(white, ar1, 'white noise should exceed AR(1)');
tc.verifyGreaterThan(ar1, sine, 'AR(1) should exceed a pure sine');
end

% ------------------------------------------------------------------
% Symbolic entropy needs a binary series. It previously failed inside
% bin2dec, whose error names neither the argument nor the requirement.
% ------------------------------------------------------------------
function testSymbolicRejectsNonBinaryInputClearly(tc)
tc.verifyError(@() ent(tc.TestData.ar1, algorithm="symbolic"), 'ent:notBinary');
tc.verifyError(@() ent_symbolic(tc.TestData.ar1, 3), 'ent_symbolic:notBinary');
end

function testSymbolicWorksOnBinaryInput(tc)
b = double(tc.TestData.ar1 > median(tc.TestData.ar1));
v = ent(b, algorithm="symbolic", dim=3);
tc.verifyGreaterThan(v, 0, 'symbolic entropy of a thresholded series should be positive');
tc.verifyLessThanOrEqual(v, 1 + 1e-9, 'the normalised measure cannot exceed 1');
end

% ------------------------------------------------------------------
% Wrapper dispatch and agreement with the direct calls.
% ------------------------------------------------------------------
function testWrapperMatchesDirectCalls(tc)
x = tc.TestData.ar1;
tc.verifyEqual(ent(x), ent_samp(x,2,0.2), 'RelTol', 1e-12);
tc.verifyEqual(ent(x, algorithm="approximate"), ent_ap(x,2,0.2), 'RelTol', 1e-12);
tc.verifyEqual(ent(x, algorithm="permutation", dim=3), ent_permu(x,3,1), 'RelTol', 1e-12);
end

function testMultiscaleReturnsOnePerScale(tc)
scale = 8;
[v, extra] = ent(tc.TestData.ar1, algorithm="multiscale", scale=scale);
tc.verifyEqual(numel(v), scale, sprintf( ...
    'multiscale returned %d values for scale=%d', numel(v), scale));
tc.verifyTrue(all(isfinite(v)), 'multiscale produced non-finite entropies');
tc.verifyEqual(extra.scale, scale);
end

function testCrossEntropyNeedsTwoSeries(tc)
x = tc.TestData.ar1; y = tc.TestData.white;
tc.verifyGreaterThan(ent(x, y=y, algorithm="sample"), 0);
tc.verifyGreaterThan(ent(x, y=y, algorithm="approximate"), 0);
tc.verifyError(@() ent(x, y=y, algorithm="permutation"), 'ent:noCrossVariant');
tc.verifyError(@() ent(x, y=y(1:10)), 'ent:lengthMismatch');
end

function testInvalidInputsAreRejected(tc)
x = tc.TestData.ar1;
tc.verifyError(@() ent(x, algorithm="wavelet"), 'ent:unknownAlgorithm');
xn = x; xn(3) = NaN;
tc.verifyError(@() ent(xn), 'ent:nanInput');
tc.verifyError(@() ent([1;2;3], dim=9), 'ent:tooShort');
end

% ------------------------------------------------------------------
% No console output. Two functions in this family printed on every call.
% ------------------------------------------------------------------
function testFamilyIsSilent(tc)
x = nonantest.signals('ar1', 400, 0.5);
b = double(x > median(x));
calls = { ...
    'ent_samp',     @() ent_samp(x, 2, 0.2); ...
    'ent_symbolic', @() ent_symbolic(b, 3); ...
    'ent_permu',    @() ent_permu(x, 3, 1); ...
    'ent',          @() ent(x)};
for k = 1:size(calls,1)
    out = evalc('calls{k,2}();');
    tc.verifyEmpty(strtrim(out), sprintf( ...
        '%s wrote to the console: "%s"', calls{k,1}, strtrim(out)));
end
end

function testRunsHeadless(tc)
s = nonantest.sideEffects(@() ent(tc.TestData.white));
tc.verifyFalse(s.errored, 'ent errored');
tc.verifyEqual(s.figures, 0, 'ent opened a figure');
end
