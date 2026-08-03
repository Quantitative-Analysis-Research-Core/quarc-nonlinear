function [ent_samp_value] = ent_samp(x, dim, radius, flag)
%SAMPEN Calculate the sample entropy of a time series.
%
%   ent_samp_value = ent_samp(x, dim, radius, flag)
%
%   Inputs:
%       x - A vector containing the time series x.
%       dim    - Embedding dimension.
%       radius    - Tolerance threshold. If flag is 'prop', then radius is a proportion
%              of the standard deviation of x (e.g., 0.2 means 0.2*std(x)).
%              If flag is 'const', then radius is used as the constant threshold.
%       flag - (optional) A string that specifies how to interpret radius.
%              Use 'prop' (default) if radius is a proportion of std(x),
%              or 'const' if radius is a constant.
%
%   Output:
%       ent_samp_value - The computed sample entropy value.
%
%   The sample entropy is defined as:
%       SampEn = -log( A / B )
%   where:
%       B = number of pairs of vectors of length dim that are similar.
%       A = number of pairs of vectors of length dim+1 that are similar.
%
%   Reference:
%       Richman, J. S. & Moorman, J. R. (2000), 
%       "Physiological time-series analysis using approximate entropy and sample entropy",
%       American Journal of Physiology-Heart and Circulatory Physiology, 278(6), H2039-H2049.
%
%   Written by Aaron D. Likens and Seung Kyeom Kim


% Ensure that the x is a column vector
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

x = x(:);
N = length(x);

% Check if the time series is long enough
if N < dim+1
    error('The time series must have at least m+1 data points.');
end

% Set default radius and flag if not provided
if nargin < 3 || isempty(radius)
    radius = 0.2;
end
if nargin < 4 || isempty(flag)
    flag = 'prop';
end

% If radius is a proportion of std(x), update radius accordingly.
if strcmpi(flag, 'prop')
    radius = radius * std(x);
end

%% Create embedding vectors of length dim
% Each row of X is a vector of dim consecutive x points.
% Richman & Moorman count both A and B over the same N-dim templates, so the
% two index the same set and -log(A/B) is a conditional probability. Using the
% N-dim+1 templates that exist for length dim would use one more template for
% B than for A and bias the ratio by (N-dim)/(N-dim+1).
X = zeros(N - dim, dim);
for i = 1:(N - dim)
    X(i, :) = x(i:i+dim-1);
end

% Count the number of similar pairs for vectors of length dim (B)
B = 0;
for i = 1:size(X, 1)
    for j = i+1:size(X, 1)
        % Using Chebyshev distance: maximum absolute difference
        if max(abs(X(i,:) - X(j,:))) <= radius
            B = B + 1;
        end
    end
end

%% Create embedding vectors of length dim+1
X1 = zeros(N - dim, dim+1);
for i = 1:(N - dim)
    X1(i, :) = x(i:i+dim);
end

% Count the number of similar pairs for vectors of length dim+1 (A)
A = 0;
for i = 1:size(X1, 1)
    for j = i+1:size(X1, 1)
        if max(abs(X1(i,:) - X1(j,:))) <= radius
            A = A + 1;
        end
    end
end

% Both counts run over the same N-dim templates, so any common normalisation
% cancels in the ratio below and is omitted.

%% Calculate sample entropy
if B == 0
    ent_samp_value = Inf;
else
    ent_samp_value = -log(A / B);
end

end