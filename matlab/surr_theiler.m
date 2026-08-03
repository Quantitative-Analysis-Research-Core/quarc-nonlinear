function z =surr_theiler(y,algorithm)
% z=surr_theiler(y,algorithm)
% inputs  - y, time series to be surrogated
%              algorithm - the type of algorithm to be completed
% outputs - z, surrogated time series
% Remarks
% - This code creates a surrogate time series according to Algorithm 0,
%   Algorithm 1 or Algorithm 2.
% Future Work
% - None.
% References
% - Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., & Doyne 
%   Farmer, J. (1992). Testing for nonlinearity in time series: the 
%   method of surrogate data. Physica D: Nonlinear Phenomena, 58(1–4), 
%   77–94. https://doi.org/10.1016/0167-2789(92)90102-S
% Jun 2015 - Modified by Ben Senderling
%          - Added function help section and plot commands for user
%            feedback
%          - The code was originally created as two algorithms. It was
%            modified so one code included both functions.
% Jul 2020 - Modified by Ben Senderling, bmchnonan@unomaha.edu
%          - Changed file and function name.
%          - Added reference.
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.
%% Begin code
switch (algorithm)
    case 0
        z=randn(size(y));
        [~,idx]=sort(z);
        z=y(idx);
    case 1
        z=surr1(y,1);
    case 2
        z=surr1(y,2);
end

end

function z = surr1(x,algorithm)

[r,c] = size(x);

y= zeros(r,c);

if abs(algorithm)==2
    ra= randn(size(x));
    [sr,~]= sort(ra);
    [sx,xi]= sort(x);
    [~,xii]= sort(xi);
    for k=1:c
        y(:,k)= sr(xii(:,k));
    end
else
    y= x;
end
m= mean(y);
y= y-m(ones(r,1),:);

fy = fft(y);

% randomizing phase
% Randomise the phase of the positive frequencies only, then mirror the
% conjugate onto the negative frequencies. This makes the rotated spectrum
% conjugate symmetric, so the inverse transform is real by construction and
% |FFT| is preserved exactly, which is the definition of a Fourier transform
% surrogate. Drawing an independent phase for every bin leaves the spectrum
% without that symmetry, so ifft returns a complex series and real() discards
% the imaginary part, losing half the energy.
% DC is left unrotated, as is the Nyquist bin when r is even, since both must
% remain real. One phase draw is shared across columns, as before, so the
% cross-spectra of a multivariate input are preserved.
nHalf = floor((r-1)/2);
rot = ones(r,1);
if nHalf > 0
    ph = 2*pi*rand(nHalf,1);
    rot(2:nHalf+1) = exp(1i*ph);
    rot(r:-1:r-nHalf+1) = conj(exp(1i*ph));
end
rot = rot(:,ones(1,c));
fyy= fy .* rot;

yy= real(ifft(fyy)) +  m(ones(r,1),:);

if abs(algorithm)==2
    [~,yyi] = sort(yy);
    [~,yyii] = sort(yyi);
    for k=1:c
        z(:,k) = sx(yyii(:,k));
    end
else
    z= yy;
end

end





