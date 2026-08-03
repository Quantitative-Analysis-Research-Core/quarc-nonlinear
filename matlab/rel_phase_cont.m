function [crp, crpH] = rel_phase_cont(data1,data2,samprate)
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

data1 = data1 - min(data1) - range(data1)/2;
data2 = data2 - min(data2) - range(data2)/2;

% Calculate a basic phase

crp = atan((diff(data1)*samprate.*data2(1:end-1) - diff(data2)*samprate.*data1(1:end-1))./(data1(1:end-1).*data2(1:end-1) - diff(data2)*samprate.*diff(data1)*samprate))*180/pi;

% Calculate the phase angle using a Hilbert tranform
env1H = hilbert(data1);
env2H = hilbert(data2);
crpH = atan((imag(env1H).*data2 - imag(env2H).*data1)./(data1.*data2 + imag(env1H).*imag(env2H)))*180/pi;

return
