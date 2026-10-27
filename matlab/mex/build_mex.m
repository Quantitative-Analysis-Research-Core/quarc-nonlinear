function build_mex()
%BUILD_MEX Compile the C++ core's MEX gateways into matlab/private.
%
%   build_mex
%
%   Each matlab/mex/quarc_*.cpp is linked with every source in cpp/src. The
%   resulting quarc_*.mex* files land in matlab/private, where only the
%   library's own functions can call them. A function with a compiled kernel
%   uses it when the MEX file exists and falls back to its MATLAB code when
%   it does not, so the library works unbuilt; the MATLAB code stays the
%   specification the kernels are tested against.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
core = fullfile(repo, 'cpp');
out = fullfile(repo, 'matlab', 'private');
src = dir(fullfile(core, 'src', '*.cpp'));
src = fullfile({src.folder}, {src.name});
gw = dir(fullfile(here, 'quarc_*.cpp'));
for g = gw'
    fprintf('building %s\n', g.name);
    mex('-R2018a', '-silent', ['-I' fullfile(core, 'include')], ...
        'CXXFLAGS=$CXXFLAGS -std=c++17 -ffp-contract=off', 'CXXOPTIMFLAGS=-O3 -DNDEBUG', ...
        '-outdir', out, fullfile(g.folder, g.name), src{:});
end
end
