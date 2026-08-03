function varargout = JRQA021525(varargin)
%JRQA021525 (deprecated) Use jrqa instead.
%   JRQA021525 has been renamed to jrqa. This shim forwards every argument and
%   output unchanged, so existing scripts keep working, and warns once per
%   session.
%
%   Silence with: warning('off','JRQA021525:deprecated')
%
%   See also JRQA.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

warning('JRQA021525:deprecated', ...
    'JRQA021525 has been renamed to jrqa. Update your code; this shim will be removed.');
[varargout{1:nargout}] = jrqa(varargin{:});
end
