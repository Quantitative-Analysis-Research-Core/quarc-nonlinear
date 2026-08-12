function s = fit_linearity(xv, yv)
%FIT_LINEARITY Whether a fitted stretch is actually straight.
%
%   s = quarctest.fit_linearity(xv, yv)
%
%   Returns a struct of diagnostics for the straight-line fit of yv on xv:
%
%     n          points used
%     r2         R^2 of the linear fit
%     curvature  quadratic coefficient, scaled by the window's own extent:
%                b2 * range(x)^2 / range(y). Dimensionless, signed, and
%                readable as the fraction of the window's rise contributed by
%                the quadratic term. Zero is straight.
%     maxDev     largest absolute residual as a fraction of range(y)
%     runsZ      z-score of the number of sign runs in the residuals against
%                what independent noise would give. Strongly negative means
%                the residuals arrive in long same-sign blocks, which is what
%                a systematic bend looks like.
%     slopeCv    standard deviation of the point-to-point local slopes,
%                divided by the fitted slope
%
%   R^2 DOES NOT MEASURE LINEARITY, WHICH IS WHY THIS EXISTS. R^2 is the
%   fraction of variance a straight line explains, and a smooth monotone bend
%   explains almost all of it. Measured on this catalogue: the Thomas system's
%   correlation-sum fit scores R^2 = 0.988 over a region that is visibly
%   curved, and returns D2 = 2.48 against a published 1.84. Rosenstein's fit
%   on the same system scores 0.977 while running well past the knee into
%   saturation. Any of these four other diagnostics separates those cases from
%   a genuinely straight one; R^2 does not.
%
%   The curvature and runs statistics are the discriminating pair. Curvature
%   states how bent the window is in units of its own rise, so it is
%   comparable across systems whose curves span wildly different ranges. The
%   runs statistic detects the same thing without assuming the departure is
%   quadratic: a straight fit through a curve produces residuals that are
%   positive in the middle and negative at both ends, giving far fewer runs
%   than chance.
%
%   NONE OF THIS MAKES AN ESTIMATE CORRECT. A perfectly straight region can
%   still yield the wrong exponent: the linear congruential generator fits at
%   R^2 = 1.0000 with negligible curvature and returns a correlation dimension
%   of 1.95 against a published value of exactly 1, because its orbit fills
%   the reconstruction rather than lying on an attractor. Linearity is
%   necessary for a scaling exponent to mean anything, not sufficient.
%
%   See also QUARCTEST.METRIC_BATTERY.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    xv (:,1) double
    yv (:,1) double
end

s = struct('n', 0, 'r2', NaN, 'curvature', NaN, 'maxDev', NaN, ...
           'runsZ', NaN, 'slopeCv', NaN);

keep = isfinite(xv) & isfinite(yv);
xv = xv(keep); yv = yv(keep);
n = numel(xv);
s.n = n;
if n < 4, return, end

dx = max(xv) - min(xv);
dy = max(yv) - min(yv);

p1 = polyfit(xv, yv, 1);
r  = yv - polyval(p1, xv);

ssTot = sum((yv - mean(yv)).^2);
if ssTot > 0
    s.r2 = 1 - sum(r.^2)/ssTot;
end
if dy > 0
    s.maxDev = max(abs(r))/dy;
end

if n >= 5 && dy > 0 && dx > 0
    p2 = polyfit(xv, yv, 2);
    s.curvature = p2(1) * dx^2 / dy;
end

% Wald-Wolfowitz runs statistic on the residual signs.
sg = sign(r); sg(sg == 0) = 1;
np = sum(sg > 0); nm = sum(sg < 0);
if np > 0 && nm > 0
    runs = 1 + sum(diff(sg) ~= 0);
    mu = 2*np*nm/n + 1;
    vr = 2*np*nm*(2*np*nm - n) / (n^2 * (n - 1));
    if vr > 0
        s.runsZ = (runs - mu)/sqrt(vr);
    end
end

ls = diff(yv)./diff(xv);
ls = ls(isfinite(ls));
if ~isempty(ls) && abs(p1(1)) > 0
    s.slopeCv = std(ls)/abs(p1(1));
end
end
