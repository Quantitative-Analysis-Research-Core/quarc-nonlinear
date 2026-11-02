function ami_validate_bins(b)
%AMI_VALIDATE_BINS Argument validator for the bins option of AMI functions.
%   Accepts empty, a positive integer bin count, or a vector of at least two
%   strictly increasing, finite bin edges.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

if isempty(b)
    return
end
if isscalar(b)
    if ~(b >= 1 && b == round(b))
        error('ami:badBins', 'bins must be a positive integer count or a vector of edges.');
    end
elseif ~isvector(b) || any(~isfinite(b)) || any(diff(b) <= 0)
    error('ami:badBins', 'bin edges must be a vector of finite, strictly increasing values.');
end
end
