function tests = testRelPhaseDisc
%TESTRELPHASEDISC rel_phase_disc's event bookkeeping on hand-worked events.
%
%   x events every 20 samples from 10, y events every 10 from 5, so y's
%   first event comes before x's. Each x cycle [10,30), [30,50), ... holds
%   two y events, 15 and 5 samples before the cycle's end, so DRP = t/T*360
%   alternates 540 and 180 degrees with n = 1, 0. The y event at 5 (before
%   the first x event) and the last one (no next y event) have no row.
tests = functiontests(localfunctions);
end

function setup(tc)
tc.TestData.vis = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
end

function teardown(tc)
close all
set(groot, 'DefaultFigureVisible', tc.TestData.vis);
end

function testYLeadsX(tc)
S1 = 10:20:210;
S2 = 5:10:215;
x = sin(2*pi*(1:220)/20); y = x;
[~, ~, D] = evalc_rpd(x, y, S1, S2);
tc.verifyEqual(height(D), 20);
tc.verifyEqual(D.t', repmat([15 5], 1, 10));
tc.verifyEqual(D.T', repmat(10, 1, 20));
tc.verifyEqual(D.n', repmat([1 0], 1, 10));
tc.verifyEqual(D.DRP', repmat([540 180], 1, 10));
end

function testSameAnswerWhetherYLeadsOrNot(tc)
% Dropping y's leading event must not change any row.
S1 = 10:20:210;
x = sin(2*pi*(1:220)/20); y = x;
[~, ~, D1] = evalc_rpd(x, y, S1, 5:10:215);
[~, ~, D2] = evalc_rpd(x, y, S1, 15:10:215);
tc.verifyEqual(D1, D2);
end

function testEmptyXCycle(tc)
% An x cycle with no y event in it adds no rows, and the y events after it
% are still measured to their own next x event (never a negative time).
S1 = [10 30 50 70 90 110 130];
S2 = [15 25 55 65 75 85 95 105 115 125 135];
x = sin(2*pi*(1:140)/20); y = x;
[~, ~, D] = evalc_rpd(x, y, S1, S2);
tc.verifyTrue(all(D.t > 0));
tc.verifyEqual(D.t', [15 5 15 5 15 5 15 5 15 5]);
tc.verifyEqual(D.n', [1 0 1 0 1 0 1 0 1 0]);
end

function [FR, PC, D] = evalc_rpd(x, y, S1, S2)
[~, FR, PC, D] = evalc('rel_phase_disc(x, y, S1, S2)');
end
