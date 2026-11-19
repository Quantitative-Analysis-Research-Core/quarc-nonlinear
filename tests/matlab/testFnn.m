function tests = testFnn
%TESTFNN fnn against a brute-force false-nearest-neighbour count.
%
%   fnn finds neighbours with a k-d tree of 512-point leaves. Until 2026-10
%   its ball-within-bounds test was inverted, ending the search as soon as
%   the ball reached a cell bound, so on series long enough to need more
%   than one leaf 2-10% of "nearest" neighbours were not the nearest and dE
%   was inflated. A brute-force count over all pairs has no tree to get
%   wrong, so the two must agree exactly.
tests = functiontests(localfunctions);
end

function testMatchesBruteForceAcrossLeaves(tc)
% N > 512 so the tree has several leaves; dims 1-4 so the search crosses
% cell bounds in more than one coordinate.
cases = {quarctest.signals('lorenz', 1500), 8; quarctest.signals('white', 1200), 1};
for k = 1:size(cases, 1)
    x = cases{k,1}; delay = cases{k,2};
    [~, dE] = fnn(x, delay, 4, 15, 2, 0);
    tc.verifyEqual(dE, localBruteForce(x, delay, 4, 15, 2), 'AbsTol', 1e-12, ...
        sprintf('case %d: fnn dE differs from the brute-force count', k));
end
end

function testShortSeriesUnchanged(tc)
% Below 512 points there is one leaf and the stopping rule never mattered.
x = quarctest.signals('henon', 400);
[~, dE] = fnn(x, 1, 4, 15, 2, 0);
tc.verifyEqual(dE, localBruteForce(x, 1, 4, 15, 2), 'AbsTol', 1e-12);
end

function testThresholdPicksFirstDimensionBelowIt(tc)
% thresh replaces the default rules: dim is the first dimension at or below
% it, the same with speed 0 or 1, and maxdim if none reaches it.
x = quarctest.signals('lorenz', 1500);
[~, dE] = fnn(x, 8, 6, 15, 2, 0);
for t = [0.5 0.1 0.01]
    want = find(dE <= t, 1);
    if isempty(want), want = 6; end
    [~, d0] = evalc('fnn(x, 8, 6, 15, 2, 0, t)');
    [~, d1] = evalc('fnn(x, 8, 6, 15, 2, 1, t)');
    tc.verifyEqual(d0, want, sprintf('thresh %g, speed 0', t));
    tc.verifyEqual(d1, want, sprintf('thresh %g, speed 1', t));
end
want = find(dE <= 0, 1);
if isempty(want), want = 6; end
[~, d] = evalc('fnn(x, 8, 6, 15, 2, 0, 0)');
tc.verifyEqual(d, want, 'thresh 0');
tc.verifyError(@() fnn(x, 8, 6, 15, 2, 0, -0.1), 'fnn:badThresh');
tc.verifyError(@() fnn(x, 8, 6, 15, 2, 0, 2), 'fnn:badThresh');
tc.verifyError(@() fnn(x, 8, 6, 15, 2, 0, [0.1 0.2]), 'fnn:badThresh');
end

function dE = localBruteForce(x, delay, maxdim, Rtol, Atol)
n = numel(x) - delay*maxdim;
RA = std(x);
dE = zeros(maxdim, 1);
y = [];
for j = 1:maxdim
    y = [y, x(1+delay*(j-1):n+delay*(j-1))]; %#ok<AGROW>
    z = x(1+delay*j:n+delay*j);
    D = sqrt(max(sum(y.^2, 2) + sum(y.^2, 2)' - 2*(y*y'), 0));
    D(1:n+1:end) = Inf;
    [dn, nn] = min(D, [], 2);
    % recompute the winning distance directly, as fnn does
    dn = sqrt(sum((y - y(nn,:)).^2, 2));
    dz = z - z(nn);
    L = abs(dz) > dn*Rtol | sqrt(dn.^2 + dz.^2)/RA > Atol;
    dE(j) = sum(L)/n;
end
end
