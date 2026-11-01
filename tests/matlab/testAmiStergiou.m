function tests = testAmiStergiou
%TESTAMISTERGIOU The deprecated AMI_Stergiou entry point.
%
%   ami_stergiou is kept unchanged so it reproduces published results, so
%   these tests pin what it returns rather than what it should return. Its
%   replacement, ami(x, L, Algorithm="histogram"), fixes a bin off-by-one;
%   the deprecation warning says values move by ~0.4% and tau does not.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
tc.TestData.warn = warning('off', 'AMI_Stergiou:deprecated');
tc.TestData.lorenz = quarctest.signals('lorenz', 2000);
tc.TestData.white  = quarctest.signals('white', 2000);
end

function teardownOnce(tc)
warning(tc.TestData.warn);
end

function testWarnsDeprecated(tc)
warning('on', 'AMI_Stergiou:deprecated');
c = onCleanup(@() warning('off', 'AMI_Stergiou:deprecated'));
tc.verifyWarning(@() ami_stergiou(tc.TestData.lorenz, 30), 'AMI_Stergiou:deprecated');
end

function testCurveLayout(tc)
L = 40;
[~, curve] = ami_stergiou(tc.TestData.lorenz, L);
tc.verifySize(curve, [L+1, 2]);
tc.verifyEqual(curve(:,1), (0:L)');
tc.verifyGreaterThanOrEqual(curve(:,2), 0);
tc.verifyEqual(max(curve(:,2)), curve(1,2), 'lag 0 should carry the most information');
end

function testFirstMinimumMatchesReplacement(tc)
L = 60;
[tau, curve] = ami_stergiou(tc.TestData.lorenz, L);
[tauNew, curveNew] = ami(tc.TestData.lorenz, L, Algorithm="histogram");
tc.assertNotEmpty(tau);
tc.verifyEqual(tau(1,1), tauNew);
tc.verifyEqual(curve(:,2), curveNew(:,2), 'RelTol', 0.02);
end

function testTauRowsAreLocalMinima(tc)
[tau, curve] = ami_stergiou(tc.TestData.lorenz, 60);
v = curve(:,2);
for r = 1:size(tau, 1) - 1          % the last row is the 20% fallback
    i = tau(r,1) + 1;
    tc.verifyTrue(v(i-1) >= v(i) && v(i) <= v(i+1), sprintf('row %d', r));
end
end

function testTwoSeriesFormIndependentNearZero(tc)
x = tc.TestData.white;
y = circshift(x, 997);
tc.verifyLessThan(ami_stergiou(x, y), 0.1 * ami_stergiou(x, x));
end

function testTwoSeriesFormIsSymmetric(tc)
x = tc.TestData.lorenz;
y = circshift(x, 9);
tc.verifyEqual(ami_stergiou(x, y), ami_stergiou(y, x), 'RelTol', 1e-12);
end
