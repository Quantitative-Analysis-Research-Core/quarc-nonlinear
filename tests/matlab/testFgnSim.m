function tests = testFgnSim
%TESTFGNSIM fgn_sim against the analytic autocovariance of fractional
%   Gaussian noise.
%
%   For fGn with Hurst exponent H and unit variance,
%
%     gamma(k) = (|k-1|^(2H) - 2|k|^(2H) + |k+1|^(2H)) / 2,
%
%   so the lag-1 autocovariance is 2^(2H-1) - 1. fgn_sim uses circulant
%   embedding (Davies-Harte), which is exact, so the sample statistics must
%   match within sampling error. The check here is independent of dfa and of
%   quarctest.signals('fgn'), so a fault in either cannot hide one here.
tests = functiontests(localfunctions);
end

function setup(tc)
tc.TestData.rng = rng;
rng(20261101, 'twister');
end

function teardown(tc)
rng(tc.TestData.rng);
end

function testLengthAndShape(tc)
for n = [2 3 100 1001]
    tc.verifySize(fgn_sim(n, 0.7), [1 n], sprintf('n = %d', n));
end
end

function testRealAndFinite(tc)
x = fgn_sim(512, 0.3);
tc.verifyTrue(isreal(x) && all(isfinite(x)));
end

function testAutocovarianceMatchesTheory(tc)
% Lags 1-3, pooled over 200 series of 1000. The true mean is 0, so the
% estimate is NOT centred: subtracting the sample mean biases it low by
% about var(mean) ~ n^(2H-2), which is 0.06 at H = 0.8 and would fail a
% correct generator.
for H = [0.2 0.5 0.8]
    c = zeros(200, 3);
    for k = 1:200
        x = fgn_sim(1000, H);
        for L = 1:3
            c(k,L) = mean(x(1:end-L) .* x(1+L:end));
        end
    end
    lag = 1:3;
    gamma = (abs(lag-1).^(2*H) - 2*abs(lag).^(2*H) + abs(lag+1).^(2*H)) / 2;
    tc.verifyEqual(mean(c), gamma, 'AbsTol', 0.02, sprintf('H = %.1f', H));
end
end

function testUnitVarianceZeroMean(tc)
x = cell2mat(arrayfun(@(~) fgn_sim(1000, 0.6), (1:100)', 'UniformOutput', false));
tc.verifyEqual(mean(x(:)), 0, 'AbsTol', 0.05);
tc.verifyEqual(mean(x(:).^2), 1, 'RelTol', 0.05);   % uncentred, as above
end

function testSeededOutputIsReproducible(tc)
rng(7); a = fgn_sim(300, 0.75);
rng(7); b = fgn_sim(300, 0.75);
tc.verifyEqual(a, b);
end
