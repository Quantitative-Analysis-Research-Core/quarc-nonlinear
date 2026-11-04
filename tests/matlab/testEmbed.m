function tests = testEmbed
%TESTEMBED embed (Judd's delay embedding) against hand-built matrices.
%
%   embed returns one ROW per lag, ordered as sort(lags), and one column per
%   embed time t, with row i holding z(t - lags(i)). With lags 0:lag:..., row
%   1 is therefore the most recent sample and the last row the oldest: the
%   reverse of psr's column order. These tests pin that layout, since callers
%   index into it directly.
tests = functiontests(localfunctions);
end

function testDimLagMatchesExplicitLags(tc)
z = (1:20)';
tc.verifyEqual(embed(z, 4, 3), embed(z, 0:3:9));
end

function testLayoutByHand(tc)
z = 1:7;
x = embed(z, 3, 2);
% lags [0 2 4], m = 7 - 4 = 3 columns, t = (1:3) + 4
tc.verifyEqual(x, [5 6 7; 3 4 5; 1 2 3]);
end

function testDefaultsAreDim3Lag1(tc)
z = rand(1, 30);
tc.verifyEqual(embed(z), embed(z, 3, 1));
end

function testColumnAndRowInputsAgree(tc)
z = rand(25, 1);
tc.verifyEqual(embed(z, 3, 2), embed(z', 3, 2));
end

function testSizeIsDimByNMinusWindow(tc)
n = 100; dim = 5; lag = 7;
tc.verifySize(embed(rand(n, 1), dim, lag), [dim, n - (dim-1)*lag]);
end

function testUnsortedLagsAreSorted(tc)
z = rand(1, 40);
tc.verifyEqual(embed(z, [4 0 2]), embed(z, [0 2 4]));
end

function testTwoOutputsSplitOffOneStepAhead(tc)
% With two outputs and no negative lag, a lag of -1 (one step into the
% future) is added and returned separately as y.
z = 1:10;
[x, y] = embed(z, 3, 1);
tc.verifyEqual(x, [3:9; 2:8; 1:7]);
tc.verifyEqual(y, 4:10);
end

function testNegativeLagsGoToSecondOutput(tc)
z = 1:10;
[x, y] = embed(z, [-2 0 1]);
% window 3, t = 2..8: rows z(t+2), z(t), z(t-1)
tc.verifyEqual(y, 4:10);
tc.verifyEqual(x, [2:8; 1:7]);
end

function testMatrixInputErrors(tc)
tc.verifyError(@() embed(rand(3, 3), 2, 1), ?MException);
end
