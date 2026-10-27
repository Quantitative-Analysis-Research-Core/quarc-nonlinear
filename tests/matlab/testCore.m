function tests = testCore
%TESTCORE The compiled kernels (cpp/, via matlab/mex) against the MATLAB code.
%
%   Each function with a kernel runs it when the MEX file is built and its
%   own MATLAB code otherwise; QUARC_NO_MEX=1 forces the MATLAB code. Here
%   both paths run on the same inputs, including rounded data with many
%   tied distances, where a kernel that merely found *a* nearest neighbour
%   rather than MATLAB's would differ. Results must be identical, except
%   where the answer passes through log(): the C library's log and MATLAB's
%   can differ in the last bit (sample entropy, the Rosenstein curve).
%   Skipped when the kernels are not built (matlab/mex/build_mex.m).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
names = {'ent_samp', 'line_hist', 'corr_dim', 'fnn', 'lye_r'};
built = cellfun(@(n) isfile(fullfile(localRepo(), 'matlab', 'private', ...
    ['quarc_' n '.' mexext])), names);
tc.assumeTrue(all(built), 'compiled kernels not built; run matlab/mex/build_mex');
tc.TestData.restore = getenv('QUARC_NO_MEX');
end

function teardownOnce(tc)
if isfield(tc.TestData, 'restore'), setenv('QUARC_NO_MEX', tc.TestData.restore); end
end

function testEntSamp(tc)
for x = localSeries()
    localSame(tc, @() [ent_samp(x{1}, 1, 0.2), ent_samp(x{1}, 2, 0.2), ent_samp(x{1}, 3, 0.15)], 1e-12);
end
end

function testLineHist(tc)
rng(3);
for n = [2 5 40 300]
    for th = [0 1 3]
        for ty = ["rqa" "crqa"]
            x = cumsum(randn(n, 1)); a = -abs(x - x');
            localSame(tc, @() localLines(a, 0.5, char(ty), th), 0);
        end
    end
end
y = cumsum(randn(120, 1)); x = cumsum(randn(120, 1));
localSame(tc, @() localLines({-abs(x - x'), -abs(y - y')}, 0.8, 'jrqa', 2), 0);
end

function testCorrDim(tc)
for x = localSeries()
    v = x{1}(1:min(end, 1200));
    localSame(tc, @() localCorrDim(v), 0);
end
end

function testFnn(tc)
for x = localSeries()
    localSame(tc, @() localFnn(x{1}), 0);
end
end

function testLyeR(tc)
for x = localSeries()
    v = x{1}(1:min(end, 1500));
    localSame(tc, @() localLyeR(v), 1e-12);
end
end

% ---------------------------------------------------------------- helpers

function localSame(tc, f, tol)
setenv('QUARC_NO_MEX', '1');  spec = f();
setenv('QUARC_NO_MEX', '0');  fast = f();
if tol == 0
    tc.verifyTrue(isequaln(fast, spec), 'kernel and MATLAB code differ');
else
    tc.verifyEqual(fast, spec, 'AbsTol', tol);
end
end

function out = localLines(a, r, ty, th)
[rec, d, v, A] = line_hist([], a, r, ty, th);
out = {rec, d, v, A};
end

function out = localCorrDim(x)
[c, e] = corr_dim(x, 3, 3, 0);
out = {c, e};
end

function out = localFnn(x)
[d, dE] = fnn(x, 3, 10, 15, 2, 0);
out = [d; dE];
end

function out = localLyeR(x)
[s, l, o] = lye_r(x, 100, 3, 4, [0 0.5 1 2], 0.1, 0);
out = [s; l; o(:)];
end

function c = localSeries()
rng(7);
c = {quarctest.signals('lorenz', 2048), quarctest.signals('henon', 512), ...
     quarctest.signals('white', 1500), round(quarctest.signals('lorenz', 1500)*2)/2, ...
     round(randn(1500, 1)*3)};
end

function r = localRepo()
r = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
