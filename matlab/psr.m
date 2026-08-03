function xpsr = psr(x,delay,dim)
%   Phase space reconstruction
%   Input
%      x : Time series that needs phase space reconstruction
%      delay : Optimal time delay
%      dim : Optimal embedding dimension
%   Output
%      xpsr : M x dim*DIM matrix

% Get time series size
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

N = height(x);
DIM = width(x);

% Compute length of phase space reconstructed data
M = N-(dim-1)*delay;

% Time delay embedding
xpsr=zeros(M,dim*DIM);
for i=1:dim
    xpsr(1:M,(1:DIM)+DIM*(i-1)) = x((1+(i-1)*delay):(N-(dim-i)*delay),:);
end

end