function tests = testLineHist
%TESTLINEHIST Line tabulation against recurrence plots built by hand.
%   line_hist tabulated line lengths with hist(d) and its ten default bins.
%   With nine or fewer lines on a diagonal every integer label lands in its
%   own bin and the answer is exact, which is why short-series tests never
%   caught it; with ten or more, adjacent labels merge, their lengths sum
%   into fictitious long lines, and the separate zeros-bin drop removed real
%   lines along with the zeros. On a 2000-sample Lorenz reconstruction at
%   2.5% recurrence that undercounted lines by 20% and inflated DET, MeanL
%   and MaxL. These tests pin the tabulation to recurrence plots whose line
%   lengths are known by construction, on both sides of the ten-line
%   boundary.
tests = functiontests(localfunctions);
end

function teardown(~)
dbclear all
end

% ------------------------------------------------------------------
function testFewLinesPerDiagonalIsExact(tc)
lens = [1 2 3 2 1];
tc.verifyEqual(sort(diagLens(makeRP(40, lens))), sort(lens), ...
    'line lengths wrong in the few-lines regime');
end

function testManyLinesPerDiagonalIsExact(tc)
% Twelve lines on one diagonal: the regime hist(d)'s ten bins could not
% represent. Before the fix this returned nine lines including two of
% length 5 that exist nowhere in the plot.
lens = [1 2 3 1 2 3 1 2 3 1 2 3];
tc.verifyEqual(sort(diagLens(makeRP(40, lens))), sort(lens), ...
    'line lengths wrong in the many-lines regime');
end

function testFullyRecurrentOffDiagonalIsKept(tc)
% A completely recurrent off-diagonal has no zeros to drop; the old
% zeros-bin heuristic deleted the whole line. Its length is n-1 and it must
% survive (the LOI, length n, is removed as always).
n = 12;
B = eye(n);
for i = 1:n-1
    B(i, i+1) = 1;
end
tc.verifyEqual(diagLens(B), n - 1, 'full off-diagonal line lost');
end

function testVerticalTabulationIsExact(tc)
% Same construct in the vertical pass: one column carrying twelve known
% vertical lines.
lens = [1 2 3 1 2 3 1 2 3 1 2 3];
n = 40;
B = zeros(n);
c = zeros(n, 1); pos = 1;
for L = lens
    c(pos:pos+L-1) = 1; pos = pos + L + 1;
end
B(:, 3) = c;
a = B * 2 - 2;
[~, ~, vert, ~] = line_hist([], a, 1, 'rqa');
vert = vert(vert > 0);
tc.verifyEqual(sort(vert(:)'), sort(lens), 'vertical line lengths wrong');
end

function testTheilerWindowExcludesNearDiagonalLines(tc)
% A line on diagonal k=3 is tangential shadowing when the window is 5 and a
% recurrence when the window is 2. The far line on k=8 survives both.
n = 20;
B = eye(n);
for i = 1:4,  B(i, i+3) = 1; end       % length-4 line on k=3
for i = 1:3,  B(i, i+8) = 1; end       % length-3 line on k=8
a = B * 2 - 2;
[~, h5, ~, ~] = line_hist([], a, 1, 'rqa', 5);
h5 = h5(h5 > 0);
tc.verifyEqual(sort(h5(:)'), 3, 'window of 5 should leave only the far line');
[~, h2, ~, ~] = line_hist([], a, 1, 'rqa', 2);
h2 = h2(h2 > 0);
tc.verifyEqual(sort(h2(:)'), [3 4], 'window of 2 should keep both lines');
end

function testTheilerRecurrenceCountsOutsideTheBandOnly(tc)
% One recurrent pair inside the band, one outside, on n=10: with theiler=2
% only the outside pair (and its transpose) can count, over the
% n^2 - (n + 2*(n-1) + 2*(n-2)) = 56 admissible cells.
n = 10;
B = eye(n);
B(1, 2) = 1; B(2, 1) = 1;              % inside the band
B(1, 5) = 1; B(5, 1) = 1;              % outside it
a = B * 2 - 2;
[rec, ~, ~, ~] = line_hist([], a, 1, 'rqa', 2);
tc.verifyEqual(rec, 100 * 2 / 56, 'AbsTol', 1e-12);
end

function testTheilerZeroReproducesOldBehaviour(tc)
% The default must be byte-identical to the historical output: LOI excluded,
% everything else counted, %REC over n^2 - n cells.
lens = [1 2 3 2 1];
B = makeRP(40, lens);
a = B * 2 - 2;
[rec0, h0, v0, ~] = line_hist([], a, 1, 'rqa');
[recW, hW, vW, ~] = line_hist([], a, 1, 'rqa', 0);
tc.verifyEqual({rec0, h0, v0}, {recW, hW, vW});
tc.verifyEqual(rec0, 100 * sum(lens) * 1 / (40^2 - 40), 'AbsTol', 1e-12);
end

% ------------------------------------------------------------------
function B = makeRP(n, lens)
% Identity (the LOI) plus the requested line lengths along the first
% super-diagonal, separated by single gaps.
B = eye(n);
c = zeros(n - 1, 1); pos = 1;
for L = lens
    c(pos:pos+L-1) = 1; pos = pos + L + 1;
end
assert(pos - 2 <= n - 1, 'pattern does not fit');
for i = 1:n-1
    B(i, i+1) = c(i);
end
end

function h = diagLens(B)
% line_hist takes a negative distance matrix: 0 where recurrent, -2 where
% not, with radius 1, thresholds back to exactly B.
a = B * 2 - 2;
[~, h, ~, ~] = line_hist([], a, 1, 'rqa');
h = h(h > 0);
h = h(:)';
end
