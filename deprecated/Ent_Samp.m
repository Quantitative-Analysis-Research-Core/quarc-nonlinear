function varargout = Ent_Samp(varargin)
%ENT_SAMP (deprecated) Use ent_samp instead.
%   Ent_Samp has been renamed to ent_samp. This shim forwards every argument and
%   output unchanged, so existing scripts keep working, and warns once per
%   session.
%
%   Silence with: warning('off','Ent_Samp:deprecated')
%
%   See also ENT_SAMP.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

warning('Ent_Samp:deprecated', ...
    'Ent_Samp has been renamed to ent_samp. Update your code; this shim will be removed.');
[varargout{1:nargout}] = ent_samp(varargin{:});
end
