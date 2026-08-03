function xSE = ent_xsamp(x,y,dim,radius,norm)
% xSE = ent_xsamp(x,y,dim,radius,norm)
% Inputs - x, first data series
%        - y, second data series
%        - dim, vector length for matching (usually 2 or 3)
%        - radius, radius tolerance to find matches (as a proportion of the average 
%             of the SDs of the data sets, usually between 0.15 and 0.25)
%        - norm, normalization to perform
%          - 1 = max rescale/unit interval (data ranges in value from 0 - 1
%            ) Most commonly used for RQA.
%          - 2 = mean/Zscore (used when data is more variable or has 
%            outliers) normalized data has SD = 1. This is best for cross 
%            sample entropy.
%          - Set to any value other than 1 or 2 to not normalize/rescale 
%            the data
% Remarks
% - Function to calculate cross sample entropy for 2 data series using the
%   method described by Richman and Moorman (2000).
% Sep 2015 - Created by John McCamley, unonbcf@unomaha.edu
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.
%%
% Check both sets of data are the same length
xl = length(x);
yl = length(y);
if xl ~= yl
    disp('The data series need to be the same length!')
end
N = length(x);
% normalize the data ensure data fits in the same "space"
if norm == 1 %normalize data to have a range 0 - 1
    xn = (x - min(x))/(max(x) - min(x));
    yn = (y - min(y))/(max(y) - min(y));
    r = radius * ((std(xn)+std(yn))/2);
elseif norm == 2 % normalize data to have a SD = 1, and mean = 0
    xn = (x - mean(x))/std(x);
    yn = (y - mean(y))/std(y);
    r = radius;
else disp('These data will not be normalized')
end

for i = 1:N-dim
    for k = 1:dim+1
        dij(:,k) = abs(xn(1+k-1:N-dim+k-1)-yn(i+k-1));
    end
    dj = max(dij(:,1:dim),[],2);
    dj1 = max(dij,[],2);
    d = find(dj<=r);
    d1 = find(dj1<=r);
    nm = length(d);
    Bm(i) = nm/(N-dim);
    nm1 = length(d1);
    Am(i) = nm1/(N-dim);
end

Bmr = sum(Bm)/(N-dim);
Amr = sum(Am)/(N-dim);

xSE = -log(Amr/Bmr);
end