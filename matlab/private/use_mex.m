function tf = use_mex(name)
%USE_MEX True when the compiled kernel NAME is built and not switched off.
%   Setting the environment variable QUARC_NO_MEX=1 makes every function run
%   its MATLAB code, the specification, instead; tests/matlab/testCore.m
%   compares the two paths.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

tf = exist(name, 'file') == 3 && ~strcmp(getenv('QUARC_NO_MEX'), '1');
end
