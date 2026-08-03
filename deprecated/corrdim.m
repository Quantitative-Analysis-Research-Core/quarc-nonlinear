function varargout = corrdim(varargin)
%CORRDIM (deprecated) Use corr_dim instead.
%   corrdim has been renamed to corr_dim. This shim forwards every argument and
%   output unchanged, so existing scripts keep working, and warns once per
%   session.
%
%   Silence with: warning('off','corrdim:deprecated')
%
%   See also CORR_DIM.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

warning('corrdim:deprecated', ...
    'corrdim has been renamed to corr_dim. Update your code; this shim will be removed.');
[varargout{1:nargout}] = corr_dim(varargin{:});
end
