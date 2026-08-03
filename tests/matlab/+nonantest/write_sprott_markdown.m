function write_sprott_markdown(outfile)
%WRITE_SPROTT_MARKDOWN Regenerate the Sprott Appendix A reference table.
%   nonantest.write_sprott_markdown() writes
%   tests/fixtures/sprott_appendix_a.md from nonantest.sprott_catalog, so the
%   human-readable reference and the values the tests actually use cannot
%   drift apart.
%
%   Run this after changing the catalogue. It is not run by the test suite:
%   a document the tests can rewrite is not a reference.
%
%   See also NONANTEST.SPROTT_CATALOG.

if nargin < 1
    here = fileparts(mfilename('fullpath'));
    repo = fileparts(fileparts(fileparts(here)));
    outfile = fullfile(repo, 'tests', 'fixtures', 'sprott_appendix_a.md');
end

c = nonantest.sprott_catalog();
fid = fopen(outfile, 'w');
closer = onCleanup(@() fclose(fid));

w = @(varargin) fprintf(fid, varargin{:});

w('# Sprott (2003) Appendix A, reference values\n\n');
w('Transcribed from J. C. Sprott, *Chaos and Time-Series Analysis*, Oxford\n');
w('University Press (2003), Appendix A, "Common chaotic systems".\n\n');
w('This is the human-readable record of the values encoded in\n');
w('`tests/matlab/+nonantest/sprott_catalog.m`, and is generated from that\n');
w('catalogue so the two cannot drift. Regenerate with:\n\n');
w('```\nmatlab -batch "addpath(''tests/matlab''); nonantest.write_sprott_markdown"\n```\n\n');

w('## What the appendix says about these numbers\n\n');
w('> All Lyapunov exponents are base-e and were calculated using the methods\n');
w('> in Chapter 5. The Kaplan-Yorke dimension is given by eqn (5.29). The\n');
w('> correlation dimension uses at least 2 x 10^12 pairs, corresponding to\n');
w('> data sets of over two million points for each case, and have been\n');
w('> extrapolated to the limit of zero size scale.\n\n');
w('and: "Most of the results are original calculations, and all have been\n');
w('independently verified."\n\n');
w('Two consequences for this library:\n\n');
w('- lambda is in **nats** per unit time. `lye_w` accumulates log2 and returns\n');
w('  bits, so multiply its output by ln 2 before comparing.\n');
w('- D2 carries real uncertainty and is not a tight target. Values shown\n');
w('  without a tolerance are given as exact in the appendix.\n\n');

w('## Systems\n\n');
w('| section | name | category | lambda (nats) | tier | D2 | usable |\n');
w('|---|---|---|---|---|---|---|\n');
for i = 1:numel(c)
    s = c(i);
    if isfinite(s.lambda)
        lam = sprintf('%.4f', s.lambda);
    else
        lam = 'infinite';
    end
    if s.d2_err > 0
        d2 = sprintf('%.3f +/- %.3f', s.d2, s.d2_err);
    else
        d2 = sprintf('%.3f (exact)', s.d2);
    end
    w('| %s | `%s` | %s | %s | %s | %s | %s |\n', ...
        s.section, s.name, strrep(s.category, '_', ' '), lam, s.tier, d2, ...
        string(s.usable));
end

w('\n## Systems excluded from benchmarking\n\n');
for i = 1:numel(c)
    s = c(i);
    if ~s.usable
        note = strtrim(strjoin(cellstr(s.note), ' '));
        note = regexprep(note, '\s+', ' ');
        w('**`%s`** (%s)\n\n%s\n\n', s.name, s.section, note);
    end
end

w('## Counts\n\n');
cats = unique([c.category], 'stable');
w('| category | systems |\n|---|---|\n');
for k = 1:numel(cats)
    w('| %s | %d |\n', strrep(cats(k), '_', ' '), sum([c.category] == cats(k)));
end
w('| **total** | **%d** |\n\n', numel(c));
w('%d of %d are usable as benchmarks. %d have lambda in closed form, marked\n', ...
    sum([c.usable]), numel(c), sum([c.tier] == "exact"));
w('tier `exact`; the rest are well-converged numerical results.\n');

fprintf('wrote %s\n', outfile);
end
