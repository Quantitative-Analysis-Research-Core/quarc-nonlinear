function tests = testEntSamp
%TESTENTSAMP Contract and known-answer tests for ent_samp.
tests = functiontests(localfunctions);
end

function teardown(~)
dbclear all
end

% ------------------------------------------------------------------
% Definition. Richman & Moorman count A and B over the same N-dim
% templates, so the two index the same set and -log(A/B) is a conditional
% probability. Counting B over the N-dim+1 templates that exist for length
% dim uses one more template for B than for A and biases the ratio by
% (N-dim)/(N-dim+1) -- measured as +3.0e-3 at N=512, dim=2.
% ------------------------------------------------------------------
function testMatchesRichmanMoormanDefinition(tc)
x = nonantest.signals('ar1', 512, 0.7);
dim = 2; radius = 0.2 * std(x); N = numel(x);

% Reference computed directly from the definition, independently of the
% implementation under test.
X  = zeros(N-dim, dim);
X1 = zeros(N-dim, dim+1);
for i = 1:(N-dim)
    X(i,:)  = x(i:i+dim-1);
    X1(i,:) = x(i:i+dim);
end
countPairs = @(M) sum(arrayfun(@(i) ...
    sum(max(abs(M(i+1:end,:) - M(i,:)), [], 2) <= radius), 1:size(M,1)-1));
expected = -log(countPairs(X1) / countPairs(X));

tc.verifyEqual(ent_samp(x, dim, 0.2), expected, 'RelTol', 1e-12, ...
    'ent_samp must equal -log(A/B) with A and B over the same N-dim templates');
end

function testTemplateCountsAreEqual(tc)
% The defect this guards: B was built over N-dim+1 templates while A used
% N-dim, so the two counts ranged over different sets.
here = fileparts(mfilename('fullpath'));
src = fileread(fullfile(fileparts(fileparts(here)), 'matlab', 'ent_samp.m'));
src = regexprep(src, '%.*', '');
n = numel(regexp(src, 'zeros\(N - dim,', 'once')) ...
  + numel(regexp(src, 'zeros\(N - dim \+ 1,', 'once'));
tc.verifyEmpty(regexp(src, 'zeros\(N - dim \+ 1,', 'once'), ...
    ['ent_samp allocates N-dim+1 templates. A and B must range over the ' ...
     'same N-dim templates for -log(A/B) to be a conditional probability.']);
end

% ------------------------------------------------------------------
% Known-answer behaviour.
% ------------------------------------------------------------------
function testWhiteNoiseIsMoreIrregularThanAutocorrelated(tc)
% Sample entropy measures irregularity. Independent samples are maximally
% irregular; an AR(1) process is predictable from its own past.
white = ent_samp(nonantest.signals('white', 2000), 2, 0.2);
ar1   = ent_samp(nonantest.signals('ar1', 2000, 0.9), 2, 0.2);
tc.verifyGreaterThan(white, ar1, sprintf( ...
    'white noise scored %.4f, AR(1) phi=0.9 scored %.4f', white, ar1));
end

function testDecreasesWithPredictability(tc)
% Sweeping the AR coefficient sweeps predictability, so entropy should fall
% monotonically as phi rises.
phis = [0.0 0.3 0.6 0.9];
got = arrayfun(@(p) ent_samp(nonantest.signals('ar1', 2000, p), 2, 0.2), phis);
fprintf('    [known answer] ent_samp over AR(1) phi %s\n                   %s\n', ...
    mat2str(phis), mat2str(round(got,4)));
tc.verifyTrue(all(diff(got) < 0), sprintf( ...
    'entropy did not fall monotonically with phi: %s', mat2str(round(got,4))));
end

function testPeriodicSignalHasLowEntropy(tc)
% A pure sine repeats exactly, so almost every template match extends.
sine = ent_samp(nonantest.signals('sine', 2000, 25), 2, 0.2);
white = ent_samp(nonantest.signals('white', 2000), 2, 0.2);
tc.verifyLessThan(sine, 0.5 * white, sprintf( ...
    'sine scored %.4f against white noise %.4f', sine, white));
end

% ------------------------------------------------------------------
% Interface.
% ------------------------------------------------------------------
function testRadiusScalesWithStandardDeviationByDefault(tc)
% The default treats radius as a proportion of the standard deviation, so the
% result is invariant to a change of units.
x = nonantest.signals('ar1', 1000, 0.5);
tc.verifyEqual(ent_samp(10*x, 2, 0.2), ent_samp(x, 2, 0.2), 'RelTol', 1e-12, ...
    'with proportional radius, ent_samp must be scale invariant');
end

function testProducesNoConsoleOutput(tc)
% A stray fprintf reporting the internal A/B ratio shipped in this function
% and polluted every batch log that called it.
x = nonantest.signals('ar1', 500, 0.5);
out = evalc('ent_samp(x, 2, 0.2);');
tc.verifyEmpty(strtrim(out), sprintf( ...
    'ent_samp wrote to the console: "%s"', strtrim(out)));
end

function testRunsHeadless(tc)
s = nonantest.sideEffects(@() ent_samp(nonantest.signals('white', 500), 2, 0.2));
tc.verifyFalse(s.errored, 'ent_samp errored');
tc.verifyEqual(s.figures, 0, 'ent_samp opened a figure');
end
