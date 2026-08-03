function varargout = LyE_W(varargin)
%LYE_W (deprecated) Use lye_w instead.
%   LyE_W has been renamed to lye_w. This shim forwards every argument and
%   output unchanged, so existing scripts keep working, and warns once per
%   session.
%
%   Silence with: warning('off','LyE_W:deprecated')
%
%   See also LYE_W.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

warning('LyE_W:deprecated', ...
    'LyE_W has been renamed to lye_w. Update your code; this shim will be removed.');
[varargout{1:nargout}] = lye_w(varargin{:});
end
