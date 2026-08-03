function [ NCSE ] = ent_symbolic( x, dim )
% [ SymEnt ] = ent_symbolic( x, dim )
% symbolicEnt Calculates the Symbolic Entropy for given data.
% Input -   x: 1-Dimensional binary array of data
%           dim: Word length
% Output -  NCSE: Normalized Corrected Shannon Entropy
% Remarks
% - This code calculates the Symbbolic Entropy value for the provided data
%   at a given word length described by - Aziz, W., Arif, M., 2006.
%   "Complexity analysis of stride interval time series by threshold
%   dependent symbolic entropy." Eur. J. Appl. Physiol. 98: 30-40.
% Jun 2017 - Created by William Denton, unonbcf@unomaha.edu
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.
%% Begin code: (Do NOT Edit)
%% Correct orientation of array.
[r,c] = size(x);
if r > c
    x = x';
end
%% Validate. bin2dec below accepts only 0 and 1, and its own error names
% neither the argument nor the requirement.
u = unique(x(:));
if ~all(ismember(u, [0 1]))
    error('ent_symbolic:notBinary', ...
        ['x must be a binary series of 0s and 1s; found %d distinct values ' ...
         'in [%g, %g]. Threshold or symbolise the series first, for example ' ...
         'x = double(x > median(x)).'], numel(u), min(u), max(u));
end

%% Convert binary values to decimal.
words = zeros(length(x)-dim+1,1);
for i = 1:length(x)-dim+1
    words(i,1) = bin2dec(num2str(x(i:i+dim-1)));
end
%% Calculate probability.
max_words = 2^dim;
for i = 1:max_words
    P(i) = sum(words == i-1)/length(words);
    H(i) = P(i)*log2(P(i));
end
H = -sum(H(~isnan(H)));
%% Normalized Corrected Shannon Entropy
So = length(unique(words));
Sm = max_words;
CSE = H+(So-1) / (2*Sm*log(2));
CSEm = -log2(1/Sm) + (Sm-1) / (2*Sm*log(2));
NCSE = CSE/CSEm;
%% Print out Symbolic Entropy Value.
end