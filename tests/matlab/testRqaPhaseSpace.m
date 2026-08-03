function tests = testRqaPhaseSpace
%TESTRQAPHASESPACE Contract tests for phase-space input to rqa/crqa/jrqa/mdrqa.
%
%   The point of accepting a phase space directly (mirroring lyapunov, see
%   testLyapunov) is that reconstruction can be done once, inspected, and
%   shared: the estimator must treat a supplied phase space exactly as it
%   would treat the identical matrix built internally, and must not silently
%   re-embed it.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
rng(20260727, 'twister');
n = 220;
tc.TestData.x = cumsum(randn(n,1));
tc.TestData.y = cumsum(randn(n,1));
tc.TestData.z = cumsum(randn(n,1));
tc.TestData.tau = 2;
tc.TestData.dim = 3;
end

function teardown(~)
dbclear all
end

% ------------------------------------------------------------------
% Phase space as an input, matches the equivalent series input.
% ------------------------------------------------------------------
function testRqaPhaseSpaceMatchesSeriesInput(tc)
x = tc.TestData.x; tau = tc.TestData.tau; dim = tc.TestData.dim;
[rpSeries, resSeries] = rqa(x, tau, dim, "Zscore", 1);
% rqa embeds with psr(data, dim, tau) -- reproduce that exactly so the
% phase space handed in is bit-identical to what rqa would have built.
Y = psr(zscore(x), dim, tau);
[rpSpace, resSpace] = rqa(Y, tau, dim, "Zscore", 0, "PhaseSpace", true);

tc.verifyEqual(rpSpace, rpSeries, 'the recurrence plot must be unchanged');
tc.verifyEqual(resSpace.REC, resSeries.REC, 'AbsTol', 1e-10);
tc.verifyEqual(resSpace.DET, resSeries.DET, 'AbsTol', 1e-10);
tc.verifyEqual(resSpace.LAM, resSeries.LAM, 'AbsTol', 1e-10);
end

function testCrqaPhaseSpaceMatchesSeriesInput(tc)
x = tc.TestData.x; y = tc.TestData.y;
tau = tc.TestData.tau; dim = tc.TestData.dim;
[rpSeries, resSeries] = crqa([x y], tau, dim, "Zscore", 1);
Y = psr(zscore([x y]), tau, dim);
[rpSpace, resSpace] = crqa(Y, tau, dim, "Zscore", 0, "PhaseSpace", true);

tc.verifyEqual(rpSpace, rpSeries, 'the cross recurrence plot must be unchanged');
tc.verifyEqual(resSpace.REC, resSeries.REC, 'AbsTol', 1e-10);
tc.verifyEqual(resSpace.DET, resSeries.DET, 'AbsTol', 1e-10);
end

function testJrqaPhaseSpaceMatchesSeriesInput(tc)
x = tc.TestData.x; y = tc.TestData.y;
tau = tc.TestData.tau; dim = tc.TestData.dim;
[rpSeries, resSeries] = jrqa([x y], [tau tau], [dim dim], "Zscore", 1);
Y1 = psr(zscore(x), tau, dim);
Y2 = psr(zscore(y), tau, dim);
[rpSpace, resSpace] = jrqa({Y1, Y2}, [tau tau], [dim dim], "Zscore", 0, "PhaseSpace", true);

tc.verifyEqual(rpSpace, rpSeries, 'the joint recurrence plot must be unchanged');
tc.verifyEqual(resSpace.REC, resSeries.REC, 'AbsTol', 1e-10);
tc.verifyEqual(resSpace.DET, resSeries.DET, 'AbsTol', 1e-10);
end

function testMdrqaPhaseSpaceMatchesSeriesInput(tc)
x = tc.TestData.x; y = tc.TestData.y; z = tc.TestData.z;
tau = tc.TestData.tau; dim = tc.TestData.dim;
[rpSeries, resSeries] = mdrqa([x y z], tau, dim, "Zscore", 1);
% mdrqa embeds with psr(data, dim, tau), same call order as rqa.
Y = psr(zscore([x y z]), dim, tau);
[rpSpace, resSpace] = mdrqa(Y, tau, dim, "Zscore", 0, "PhaseSpace", true);

tc.verifyEqual(rpSpace, rpSeries, 'the recurrence plot must be unchanged');
tc.verifyEqual(resSpace.REC, resSeries.REC, 'AbsTol', 1e-10);
tc.verifyEqual(resSpace.DET, resSeries.DET, 'AbsTol', 1e-10);
end

% ------------------------------------------------------------------
% A supplied phase space is used verbatim, not silently re-embedded.
% ------------------------------------------------------------------
function testRqaPhaseSpaceIsNotReembedded(tc)
% If PhaseSpace were ignored and DATA re-embedded via psr(data,dim,tau),
% the M = N-(dim-1)*tau row-trimming formula would shrink the row count.
% A verbatim pass-through must keep every row.
n = 80;
Y = randn(n, 4);
[~, res] = rqa(Y, 6, 5, "Zscore", 0, "PhaseSpace", true);
tc.verifyEqual(res.Size, n, ...
    'PhaseSpace=true must use DATA verbatim, not re-embed it with TAU/DIM.');
end

function testCrqaPhaseSpaceIsNotReembedded(tc)
n = 80;
Y = randn(n, 4);
[~, res] = crqa(Y, 6, 5, "Zscore", 0, "PhaseSpace", true);
tc.verifyEqual(res.Size, n, ...
    'PhaseSpace=true must use DATA verbatim, not re-embed it with TAU/DIM.');
end

function testMdrqaPhaseSpaceIsNotReembedded(tc)
n = 80;
Y = randn(n, 4);
[~, res] = mdrqa(Y, 6, 5, "Zscore", 0, "PhaseSpace", true);
tc.verifyEqual(res.Size, n, ...
    'PhaseSpace=true must use DATA verbatim, not re-embed it with TAU/DIM.');
end

function testJrqaPhaseSpaceIsNotReembedded(tc)
n = 80;
Y1 = randn(n, 4);
Y2 = randn(n, 4);
[~, res] = jrqa({Y1, Y2}, [6 6], [5 5], "Zscore", 0, "PhaseSpace", true);
tc.verifyEqual(res.Size, n, ...
    'PhaseSpace=true must use DATA verbatim, not re-embed it with TAU/DIM.');
end

% ------------------------------------------------------------------
% Existing (non-PhaseSpace) calls are unchanged.
% ------------------------------------------------------------------
function testRqaDefaultCallStillWorksWithoutPhaseSpace(tc)
x = tc.TestData.x;
[rp, res] = rqa(x, tc.TestData.tau, tc.TestData.dim);
tc.verifyEqual(res.EMB, tc.TestData.dim);
tc.verifyEqual(res.DEL, tc.TestData.tau);
tc.verifyGreaterThan(res.REC, 0);
tc.verifyEqual(size(rp,1), res.Size);
end

function testRqaRawWideDataStillRejectedWithoutPhaseSpace(tc)
% A raw call is still required to be a single column; only PhaseSpace=true
% relaxes that shape requirement.
Y = randn(50, 3);
tc.verifyError(@() rqa(Y, 2, 3), '');
end

function testCrqaRawDataStillRequiresTwoColumns(tc)
oneCol = tc.TestData.x;
tc.verifyError(@() crqa(oneCol, 2, 3), '');
end

function testJrqaRawDataStillDoubleNotCell(tc)
x = tc.TestData.x; y = tc.TestData.y;
% Without PhaseSpace, a cell array is not an accepted raw input.
tc.verifyError(@() jrqa({x, y}, [1 1], [1 1]), 'jrqa:rawDataMustBeDouble');
end

% ------------------------------------------------------------------
% Norm="euc" cannot verify DIM against a supplied phase space -- warn
% rather than silently trust it (the pattern from lyapunov's Theiler
% window).
% ------------------------------------------------------------------
function testRqaWarnsOnEucNormWithPhaseSpace(tc)
Y = psr(zscore(tc.TestData.x), tc.TestData.dim, tc.TestData.tau);
tc.verifyWarning(@() rqa(Y, tc.TestData.tau, tc.TestData.dim, ...
    "Zscore", 0, "PhaseSpace", true, "Norm", "euc"), 'rqa:eucNormAssumesDim');
tc.verifyWarningFree(@() rqa(Y, tc.TestData.tau, tc.TestData.dim, ...
    "Zscore", 0, "PhaseSpace", true));
end

function testCrqaWarnsOnEucNormWithPhaseSpace(tc)
Y = psr(zscore([tc.TestData.x tc.TestData.y]), tc.TestData.tau, tc.TestData.dim);
tc.verifyWarning(@() crqa(Y, tc.TestData.tau, tc.TestData.dim, ...
    "Zscore", 0, "PhaseSpace", true, "Norm", "euc"), 'crqa:eucNormAssumesDim');
tc.verifyWarningFree(@() crqa(Y, tc.TestData.tau, tc.TestData.dim, ...
    "Zscore", 0, "PhaseSpace", true));
end

function testMdrqaWarnsOnEucNormWithPhaseSpace(tc)
Y = psr(zscore([tc.TestData.x tc.TestData.y tc.TestData.z]), ...
    tc.TestData.dim, tc.TestData.tau);
tc.verifyWarning(@() mdrqa(Y, tc.TestData.tau, tc.TestData.dim, ...
    "Zscore", 0, "PhaseSpace", true, "Norm", "euc"), 'mdrqa:eucNormAssumesDim');
tc.verifyWarningFree(@() mdrqa(Y, tc.TestData.tau, tc.TestData.dim, ...
    "Zscore", 0, "PhaseSpace", true));
end

function testJrqaWarnsOnEucNormWithPhaseSpace(tc)
Y1 = psr(zscore(tc.TestData.x), tc.TestData.tau, tc.TestData.dim);
Y2 = psr(zscore(tc.TestData.y), tc.TestData.tau, tc.TestData.dim);
tau2 = [tc.TestData.tau tc.TestData.tau];
dim2 = [tc.TestData.dim tc.TestData.dim];
% jrqa's Norm="euc" formula itself errors on a (1,2) DIM vector regardless
% of PhaseSpace (a pre-existing, unrelated defect); the warning must still
% fire before that error is reached.
tc.verifyWarning(@() ...
    localIgnoreLaterError(@() jrqa({Y1, Y2}, tau2, dim2, "Zscore", 0, ...
        "PhaseSpace", true, "Norm", "euc")), ...
    'jrqa:eucNormAssumesDim');
tc.verifyWarningFree(@() jrqa({Y1, Y2}, tau2, dim2, "Zscore", 0, ...
    "PhaseSpace", true));
end

function localIgnoreLaterError(fcn)
try
    fcn();
catch
    % Only the warning is under test here.
end
end

% ------------------------------------------------------------------
% Shape/type validation specific to PhaseSpace.
% ------------------------------------------------------------------
function testJrqaPhaseSpaceRequiresCellArray(tc)
Y = randn(50, 4);
tc.verifyError(@() jrqa(Y, [1 1], [1 1], "PhaseSpace", true), ...
    'jrqa:phaseSpaceMustBeCell');
end
