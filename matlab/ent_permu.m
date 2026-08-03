function [permEnt, hist] = ent_permu(x, dim, delay)
% [permEnt, hist] = ent_permu(x, dim, delay)
% inputs -  x: 1-D array of x being analyzed
%           dim: embedding dimension (order of permutation entropy) 
%           delay: time delay
% outputs - permuEnt: value calculated using a log base of 2
%           hist: number of occurences for each permutation order
% Remarks
% - It differs from the permutation entropy code found on MatLab Central in
%   one way (see MathWorks reference). The code on MatLab Central uses the 
%   log function (base e, natural log), whereas this code uses log2 (base 2
%   ), as per Bandt & Pompe, 2002. However, this code does include a lag 
%   (time delay) feature like the one on MatLab Central does.
% - Complexity parameters for time series based on comparison of 
%   neighboring values. Based on the distributions of ordinal patterns, 
%   which describe order relations between the values of a time series. 
%   Based on the algorithm described by Bandt & Pompe, 2002.
% References
% - Bandt, C., Pompe, B. Permutation entropy: A natural complexity measure 
%   for time series. Phys Rev Lett 2002, 88, 174102, 
%   doi:10.1103/PhysRevLett.88.174102
% - MathWorks: http:www.mathworks.com/matlabcentral/fileexchange/
%   37289-permutation-entropy)
% Jun 2016 - Created by Patrick Meng-Frecker, unonbcf@unomaha.edu
% Dec 2016 - Edited by Casey Wiens, email: unonbcf@unomaha.edu
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.
%%

N = length(x);  % length of time series
perm = perms(1:dim);  % create all possible permutation vectors
hist(1:length(perm)) = 0;   % designate variable to store values

for cnt1=1:N-delay*(dim-1)  % steps from 1 through length of x minus time delay multiplied by order minus 1
    [~, permVal] = sort(x(cnt1:delay:cnt1+delay*(dim-1))); % creates permutation of selected x range
    for cnt2=1:length(perm) % steps through length of possible permutation vectors
        if perm(cnt2,:) - permVal == 0  % compares current permutation of selected x to possible permutation vectors
            hist(cnt2) = hist(cnt2) + 1;    % if above comparison is equal, then adds one to bin for appropriate permutation vector
        end
    end
end

histNew = hist(hist ~= 0);  % remove any permutation orders with 0 for proper calculation
per = histNew/sum(histNew);	% ratio of each permutation vector match to total matches
permEnt = -sum(per .* log2(per));   % performs entropy calucation



