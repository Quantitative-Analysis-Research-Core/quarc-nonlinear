function varargout = linehist(varargin)
%LINEHIST (deprecated) Use line_hist instead.
%   linehist has been renamed to line_hist. This shim forwards every argument and
%   output unchanged, so existing scripts keep working, and warns once per
%   session.
%
%   Silence with: warning('off','linehist:deprecated')
%
%   See also LINE_HIST.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

warning('linehist:deprecated', ...
    'linehist has been renamed to line_hist. Update your code; this shim will be removed.');
[varargout{1:nargout}] = line_hist(varargin{:});
end
