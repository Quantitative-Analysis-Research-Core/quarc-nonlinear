function [tau, info] = ami_first_minimum(curve, fraction, caller)
%AMI_FIRST_MINIMUM Lag of the first strict local minimum of an AMI curve.
%   [TAU,INFO] = AMI_FIRST_MINIMUM(CURVE,FRACTION,CALLER) is shared by every
%   AMI estimator so that they all pick tau by the same rule. CURVE is
%   [lag, ami]. CALLER names the warning identifier, CALLER:noMinimum.
%
%   Strict on the left, so a plateau contributes at most its first point.
%   When the curve has no local minimum, tau falls back to the first lag at
%   which AMI drops to FRACTION of its lag-0 value (Abarbanel et al. 1993).

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

v = curve(:, 2);
isMin = [false; v(2:end-1) < v(1:end-2) & v(2:end-1) <= v(3:end); false];
idx = find(isMin);

info = struct();
info.allMinima = curve(idx, :);
info.usedFallback = false;
info.fractionLag = NaN;

below = find(v <= fraction * v(1), 1, 'first');
if ~isempty(below)
    info.fractionLag = curve(below, 1);
end

if ~isempty(idx)
    tau = curve(idx(1), 1);
else
    tau = info.fractionLag;              % NaN if that also failed
    info.usedFallback = true;
    if isnan(tau)
        warning([char(caller) ':noMinimum'], ...
                ['No local minimum in AMI over lags 0:%d and the curve never fell to ' ...
                 '%.0f%% of its lag-0 value. Increase maxlag, or inspect the curve directly.'], ...
                curve(end, 1), 100 * fraction);
    end
end
end
