function ref = reference_spectra(path)
%REFERENCE_SPECTRA The recomputed reference corpus, keyed by catalogue and name.
%
%   ref = quarctest.reference_spectra()
%   ref = quarctest.reference_spectra(path_to_reference_spectra_json)
%
%   Returns a containers.Map from "<source>/<lowercase name>" -- for example
%   "sprott/lorenz" or "dysts/lorenz" -- to the corpus entry for that system:
%   the full spectrum, lambdaMax, the initial condition and control parameters
%   it was computed from, the three internal checks (traceError, zeroExponent,
%   sem), and the published value it replaces.
%
%   WHY THIS EXISTS. The two catalogues the battery is scored against carried
%   exponents from two different procedures -- Sprott's trajectory separation
%   and dysts' continuous QR -- and on the 13 systems both hold, those disagree
%   by a median of 22%. python/build_reference_corpus.py recomputes every
%   system under one method; this reads the result so the MATLAB side can
%   score against it. See the corpus file's own provenance fields and
%   python/lyap_spectrum.py for the method.
%
%   Names are lower-cased on both sides of the join: the Sprott catalogue and
%   the corpus use snake_case, the dysts manifest and the corpus use CamelCase,
%   and neither should have to know which the other used.
%
%   See also QUARCTEST.RECOMPUTED_REFERENCE, QUARCTEST.SPROTT_CATALOG,
%   QUARCTEST.DYSTS_CATALOG.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    path (1,1) string = ""
end

if path == ""
    here = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % tests/
    path = fullfile(here, 'reports', 'reference_spectra.json');
end
if ~isfile(path)
    error('quarctest:referenceSpectra:noCorpus', ...
        ['No reference corpus at %s. Build one with\n' ...
         '    python/build_reference_corpus.py --out tests/reports/reference_spectra.json'], ...
        path);
end

j = jsondecode(fileread(path));
ref = containers.Map('KeyType', 'char', 'ValueType', 'any');
systems = j.systems;
% jsondecode returns a struct array when every entry has the same fields and a
% cell array otherwise; the corpus has optional fields (params, published) so
% it can come back either way.
if ~iscell(systems), systems = num2cell(systems); end
for i = 1:numel(systems)
    e = systems{i};
    key = lower(string(e.source)) + "/" + lower(string(e.name));
    ref(char(key)) = e;
end
end
