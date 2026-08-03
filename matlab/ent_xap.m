function [xAP]=ent_xap(x,y,dim,radius,k)

%[xAP]=ent_xap(x,y,dim,radius,k)
%
% inputs:    x - first time series
%            y - second time series
%            dim - something vector length
%            radius - R tolerance to find matches, proportion of the stdev
%            k - something lag
% outputs:   xAP - cross approximate entropy
%
% Remarks
% - This code finds the cross approximate entropy between two signals of
%   equal length.
%
% Future Work
% - This code should be looked over.
% - The scaling of the radius to the standard deviation may need to be
%   calculated from the average stdev of both signals and not just one.
% - The first for loop with m=dim:k:dim+k looks suspicious.
%
% Mar 2016 - Modified by Ben Senderling, email: bensenderling@gmail.com
%          - Moved the data normalization from the code that called this
%            one into this code.
%          - Changed the input radius value from a percentage to a decimal for
%            consistency with other entropy code.
%
%% Begin Code

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

x=(x-mean(x))/std(x);
y=(y-mean(y))/std(y);

N=length(x);
Cm=[];
radius=std(x)*radius;
for m=dim:k:dim+k 
	C=[];
	for i=1:(N-m+1)
		V=[x(i:m+i-1)];
		count=0;
		for j=1:(N-m+1)
			Z=[y(j:m+j-1)];
			dif=(abs(V-Z)<radius);	%two subsequences are similar if the difference between any pair 
                 				%of corrisponding measurements in the pattern is less than radius
			A=all(dif);
			count=count+A;
		end
	C=[C count/(N-m+1)];	%vector containing the 
                    		%Cim=(number of patterns similar to the one beginning at interval i)/total number  
                   		%of pattern with the same length dim
                        %display(num2str(i))
	end
	Cm=[Cm sum(C)/(N-m+1)];%vector containing the means of the Cim for subsequences of length dim and of length dim+k%
end
xAP=log(Cm(1)/Cm(2));