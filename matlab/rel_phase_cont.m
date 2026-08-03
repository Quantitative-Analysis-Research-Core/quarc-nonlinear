function [crp, crpH] = rel_phase_cont(x,y,fs)
% RELPHASE - calculates the relative phase between two segments, joint
% angles, coordinates, etc.  It is calculated using the atan2 function in
% order to preserve the quadrant in the phase portrait.  The sign of the
% relative phase is determined from segment A - segment B.  If the sign is
% negative, B leads A.  If the sign is positive, A leads B.
%
% Syntax:
%       Calculate relative phase from position and velocity data
%       rp = relphase(pa,va,pb,vb)
%           pa - angle (or position) of segment A
%           va - velocity of segment A
%           pb - angle (or position) of segment B
%           vb - velocity of segment B
%           rp - relative phase in degrees
%
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.
%% Normalize the data

x = x - min(x) - range(x)/2;
y = y - min(y) - range(y)/2;

% Calculate a basic phase

crp = atan((diff(x)*fs.*y(1:end-1) - diff(y)*fs.*x(1:end-1))./(x(1:end-1).*y(1:end-1) - diff(y)*fs.*diff(x)*fs))*180/pi;

% Calculate the phase angle using a Hilbert tranform
env1H = hilbert(x);
env2H = hilbert(y);
crpH = atan((imag(env1H).*y - imag(env2H).*x)./(x.*y + imag(env1H).*imag(env2H)))*180/pi;

return
