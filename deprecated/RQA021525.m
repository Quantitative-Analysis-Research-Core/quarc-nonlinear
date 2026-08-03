function varargout = RQA021525(varargin)
%RQA021525 (deprecated) Use rqa instead.
%   RQA021525 has been renamed to rqa. This shim forwards every argument and
%   output unchanged, so existing scripts keep working, and warns once per
%   session.
%
%   Silence with: warning('off','RQA021525:deprecated')
%
%   See also RQA.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

warning('RQA021525:deprecated', ...
    'RQA021525 has been renamed to rqa. Update your code; this shim will be removed.');
[varargout{1:nargout}] = rqa(varargin{:});
end
