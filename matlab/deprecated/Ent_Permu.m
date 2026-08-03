function varargout = Ent_Permu(varargin)
%ENT_PERMU (deprecated) Use ent_permu instead.
%   Ent_Permu has been renamed to ent_permu. This shim forwards every argument and
%   output unchanged, so existing scripts keep working, and warns once per
%   session.
%
%   Silence with: warning('off','Ent_Permu:deprecated')
%
%   See also ENT_PERMU.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

warning('Ent_Permu:deprecated', ...
    'Ent_Permu has been renamed to ent_permu. Update your code; this shim will be removed.');
[varargout{1:nargout}] = ent_permu(varargin{:});
end
