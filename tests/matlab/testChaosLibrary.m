function tests = testChaosLibrary
%TESTCHAOSLIBRARY Contract and known-answer tests for chaos_library.
%
%   The two discrete maps were the broken pair. Logistic accumulated its
%   iterate in a local x and never assigned the declared output y, so every
%   call ended in MATLAB:unassignedOutputs. Hennon used the output variable y
%   as the map's second coordinate, returning that coordinate alone as a row
%   with x discarded, and rejected the vector t its own header documents.
%
%   The known answers here are written from scratch: the logistic map at r=4
%   has the closed form x_n = sin^2(2^n asin sqrt(x0)), and henon_512.csv was
%   generated in NumPy by make_fixtures.py, independently of this library.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
tc.TestData.fx = fullfile(fileparts(fileparts(here)), 'tests', 'fixtures');
end

function teardown(~)
dbclear all
end

% ------------------------------------------------------------------
% Both maps must return something at all, and return it in the shape the
% header promises: column oriented, one t per row of y.
% ------------------------------------------------------------------
function testMapsReturnColumnOrientedOutput(tc)
cases = { ...
    'Logistic', 0.2,       3.7,       1
    'Hennon',   [0.1 0.1], [1.4 0.3], 2};
for k = 1:size(cases,1)
    name = cases{k,1};
    s = nonantest.sideEffects(@() chaos_library(name, 500, cases{k,2}, cases{k,3}));
    tc.verifyFalse(s.errored, sprintf('chaos_library(''%s'', ...) errored: %s', ...
        name, localMsg(s)));
    tc.verifyEqual(s.figures, 0, sprintf('%s opened a figure', name));

    [t, y] = chaos_library(name, 500, cases{k,2}, cases{k,3});
    tc.verifySize(y, [500 cases{k,4}], sprintf( ...
        '%s must return 500 rows and %d columns', name, cases{k,4}));
    tc.verifySize(t, [500 1], sprintf('%s must return one t per row of y', name));
    tc.verifyTrue(all(isfinite(y(:))), sprintf('%s produced non-finite values', name));
end
end

% ------------------------------------------------------------------
% Known answer: the logistic map at r = 4 is conjugate to the doubling map,
% so x_n = sin^2(2^n asin sqrt(x0)) exactly. The comparison runs 20 iterates.
% The horizon is short because the map's Lyapunov exponent of ln 2 doubles
% the rounding difference between the two routes every step: measured, the
% gap is 2e-14 at 10 iterates, 4e-11 at 20 and 1.5e-9 at 25. That growth is
% the map working, not the implementation failing.
% ------------------------------------------------------------------
function testLogisticMatchesClosedForm(tc)
x0 = 0.2;
[~, y] = chaos_library('Logistic', 20, x0, 4);
n = (0:19)';
closed = sin(2.^n * asin(sqrt(x0))).^2;
tc.verifyEqual(y, closed, 'AbsTol', 1e-9, ...
    'the logistic map at r=4 does not follow its closed form');
end

function testLogisticFixedPointsAndPeriodTwo(tc)
% Below r = 1 the origin attracts; between 1 and 3 the non-trivial fixed
% point 1-1/r attracts; just past 3 the orbit is period two.
[~, decay] = chaos_library('Logistic', 200, 0.6, 0.8);
tc.verifyLessThan(abs(decay(end)), 1e-6, 'r=0.8 must collapse to zero');

r = 2.5;
[~, fixed] = chaos_library('Logistic', 200, 0.2, r);
tc.verifyEqual(fixed(end), 1 - 1/r, 'AbsTol', 1e-9, ...
    sprintf('r=%.1f must settle on the fixed point 1-1/r', r));

[~, cyc] = chaos_library('Logistic', 400, 0.2, 3.2);
tc.verifyEqual(cyc(end), cyc(end-2), 'AbsTol', 1e-9, ...
    'r=3.2 must be period two');
tc.verifyGreaterThan(abs(cyc(end) - cyc(end-1)), 1e-3, ...
    'r=3.2 must not be a fixed point');
end

% ------------------------------------------------------------------
% Known answer: henon_512.csv is the x coordinate of the canonical map at
% a=1.4, b=0.3, written in NumPy by make_fixtures.py. Matching it tests the
% recursion, the parameter order and the column order of y at once.
%
% The run is seeded FROM the fixture rather than from the fixture's own
% initial condition. The Hennon map has a positive Lyapunov exponent near
% 0.42 per iterate, so over the fixture's 1000 discarded transients MATLAB
% and NumPy rounding separate by e^420 -- the two orbits then share an
% attractor and nothing else. Measured, no alignment of the two runs agrees
% better than 2.4, which says nothing about either implementation. Thirty
% iterates from a common state is the honest horizon. The state at fixture
% step k is (x, y) = (ref(k), b*ref(k-1)), which the map's own second
% equation supplies.
% ------------------------------------------------------------------
function testHennonReproducesTheNumpyFixture(tc)
ref = readmatrix(fullfile(tc.TestData.fx, 'henon_512.csv'));
tc.assertNumElements(ref, 512, 'fixture henon_512.csv is not 512 points');

k = 30;
[~, y] = chaos_library('Hennon', k, [ref(2) 0.3*ref(1)], [1.4 0.3]);
tc.verifyEqual(y(:,1), ref(2:k+1), 'AbsTol', 1e-8, ...
    'the Hennon x coordinate does not match the independent NumPy fixture');

% Second column is the b*x(n) coordinate, which the old code returned in
% place of the first and which no caller could reach.
tc.verifyEqual(y(2:end,2), 0.3*y(1:end-1,1), 'AbsTol', 1e-12, ...
    'the second column of y must be b*x(n)');
end

function testHennonStaysOnItsAttractor(tc)
[~, y] = chaos_library('Hennon', 2000, [0.1 0.1], [1.4 0.3]);
x = y(501:end, 1);
tc.verifyLessThan(max(abs(x)), 2, 'the Hennon orbit escaped its attractor');
tc.verifyGreaterThan(std(x), 0.1, 'the Hennon orbit collapsed to a point');
end

% ------------------------------------------------------------------
% The header documents t in vector form. For the maps that expression
% reached a colon and raised MATLAB:colon:operandsNotRealScalar.
% ------------------------------------------------------------------
function testMapsAcceptAVectorT(tc)
tvec = (0:0.5:99.5)';
for c = {{'Logistic', 0.2, 3.7, 1}, {'Hennon', [0.1 0.1], [1.4 0.3], 2}}
    a = c{1};
    s = nonantest.sideEffects(@() chaos_library(a{1}, tvec, a{2}, a{3}));
    tc.verifyFalse(s.errored, sprintf('%s rejected a vector t: %s', a{1}, localMsg(s)));
    if s.errored, continue, end

    [t, y] = chaos_library(a{1}, tvec, a{2}, a{3});
    tc.verifyEqual(t, tvec, sprintf('%s must return the sample times it was given', a{1}));
    tc.verifySize(y, [numel(tvec) a{4}], sprintf( ...
        '%s must return one row per element of t', a{1}));

    % A count and a vector of the same length describe the same run.
    [~, byCount] = chaos_library(a{1}, numel(tvec), a{2}, a{3});
    tc.verifyEqual(y, byCount, sprintf( ...
        '%s gave different orbits for a count and a vector of that length', a{1}));
end
end

% ------------------------------------------------------------------
% An unrecognised name used to fall out of the switch with both outputs
% unassigned, which reads as a fault in the solver.
% ------------------------------------------------------------------
function testUnknownSystemErrorsByName(tc)
tc.verifyError(@() chaos_library('Henon', 100, [0.1 0.1], [1.4 0.3]), ...
    'chaos_library:unknownSystem', ...
    'a misspelled system name must say so');
end

% ------------------------------------------------------------------
% A flow, to hold the shared output contract across both kinds of system.
% ------------------------------------------------------------------
function testLorenzIsUnchanged(tc)
s = nonantest.sideEffects(@() chaos_library('Lorenz', 0:0.01:20, [0 -0.01 9], [10 28 8/3]));
tc.verifyFalse(s.errored, sprintf('Lorenz errored: %s', localMsg(s)));
tc.verifyEqual(s.figures, 0, 'Lorenz opened a figure');

[t, y] = chaos_library('Lorenz', 0:0.01:20, [0 -0.01 9], [10 28 8/3]);
tc.verifySize(y, [numel(0:0.01:20) 3], 'Lorenz must return three columns');
tc.verifyEqual(numel(t), size(y,1), 'Lorenz t and y lengths disagree');
tc.verifyGreaterThan(std(y(:,1)), 1, 'the Lorenz orbit did not develop');
end

function m = localMsg(s)
if isempty(s.err), m = '(no exception)'; else, m = s.err.message; end
end
