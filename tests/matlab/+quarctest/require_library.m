function require_library()
%REQUIRE_LIBRARY Put <repo>/matlab on the path, or error trying.
%   metric_battery wraps every estimator call in a catch that records the
%   error as a fail: status on that metric. That is right for a series a
%   metric cannot handle, but when the library itself is absent it turns
%   every metric into NaN and the run completes without a word. The
%   characterize entry points call this first so a missing library fails
%   loudly at t=0 instead of silently after the whole battery.

here = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % tests/
addpath(fullfile(fileparts(here), 'matlab'));
if exist('lyapunov', 'file') ~= 2
    error('quarctest:requireLibrary:notFound', ...
          'the library at <repo>/matlab could not be put on the path');
end
end
