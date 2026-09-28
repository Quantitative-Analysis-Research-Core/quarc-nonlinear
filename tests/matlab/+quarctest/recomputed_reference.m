function c = recomputed_reference(c, source, ref)
%RECOMPUTED_REFERENCE Score a catalogue against the one-method corpus.
%
%   c = quarctest.recomputed_reference(c, "sprott")
%   c = quarctest.recomputed_reference(c, "dysts")
%   c = quarctest.recomputed_reference(c, source, ref)   % ref from reference_spectra()
%
%   Takes a catalogue struct array from quarctest.sprott_catalog or
%   quarctest.dysts_catalog and returns it with each entry's lambda replaced by
%   the recomputed value from tests/reports/reference_spectra.json. The value
%   it replaces is kept as lambdaPublished, so nothing is lost, and three
%   fields are added:
%
%     lambdaSource  where the reference now comes from, written into the
%                   summary's referenceSource column so a table can never be
%                   read without knowing which procedure scored it
%     icDependent   true where the exponent belongs to the initial condition
%                   rather than the system (conservative flows and maps with
%                   mixed phase space). The value is still supplied; a table
%                   should mark it. See sprott_systems.ic_dependent in the
%                   Python port for the reasoning.
%     lambdaChecks  the corpus entry's internal checks, so a caller can decide
%                   how much to trust the reference: traceError, zeroExponent
%                   and sem, none of which consult a published number.
%
%   A system that is not in the corpus keeps its published lambda and gets
%   lambdaSource = "published (not in corpus)", so the substitution is never
%   silent in either direction.
%
%   WHY A SEPARATE STEP RATHER THAN A CATALOGUE DEFAULT. The published values
%   are what every committed table was scored against. Swapping them inside
%   the catalogue would change every ratio in the repository with nothing in
%   the output saying so. Applying the corpus explicitly, and writing the
%   source into the row, keeps the two scorings distinguishable in the data.
%
%   See also QUARCTEST.REFERENCE_SPECTRA, QUARCTEST.CHARACTERIZE,
%   QUARCTEST.CHARACTERIZE_SERIES.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    c struct
    source (1,1) string {mustBeMember(source, ["sprott", "dysts"])}
    ref = []
end

if isempty(ref)
    ref = quarctest.reference_spectra();
end

for i = 1:numel(c)
    c(i).lambdaPublished = c(i).lambda;
    key = char(lower(source) + "/" + lower(string(c(i).name)));
    if ~isKey(ref, key)
        c(i).lambdaSource = "published (not in corpus)";
        c(i).icDependent = false;
        c(i).lambdaChecks = struct('traceError', NaN, 'zeroExponent', NaN, 'sem', NaN);
        continue
    end
    e = ref(key);
    c(i).lambda = double(e.lambdaMax);
    c(i).lambdaSource = "recomputed, one method (reference_spectra.json)";
    c(i).icDependent = isfield(e, 'icDependent') && logical(e.icDependent);
    sem = NaN;
    if isfield(e, 'sem') && ~isempty(e.sem), sem = double(e.sem(1)); end
    zero = NaN;
    if isfield(e, 'zeroExponent') && ~isempty(e.zeroExponent)
        zero = double(e.zeroExponent);
    end
    c(i).lambdaChecks = struct( ...
        'traceError',   double(e.traceError), ...
        'zeroExponent', zero, ...
        'sem',          sem);
end
end
