function [w_ent] = ent_weighted(wrp)
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

N = length(wrp); % get size of the weighted recurrence plot
for j = 1:N
    si(j) = sum(wrp(:,j)); % compute vertical weights sums
end

% Compute distribution of vertical weights sums
si_min = min(si);
si_max = max(si);

% All columns carrying equal weight gives a zero bin size, and si_min:0:si_max
% is an empty range in MATLAB, so the histogram loop below would not execute
% and p1 would never be created. A single occupied bin carries no information,
% so the entropy is zero. This arises on binary and categorical data.
if si_max == si_min
    w_ent = 0;
    return
end

bin_size = (si_max - si_min)/49; % compute bin size
S = sum(si);
% Bins are [si_min + (k-1)*bin_size, si_min + k*bin_size) for k = 1..49,
% and bin 50 holds the maximum alone: the layout the original loop over
% si_min:bin_size:si_max intended. That loop built its starts by repeated
% addition, so whether the maximum got its own bin, shared bin 49, or (when
% the last start rounded above it) fell in no bin at all depended on the
% last bit. Two distinct column sums once scored above ln 2, which is
% impossible. Assigning by index, with the maximum placed explicitly, gives
% the same bins in any language.
bin = min(floor((si - si_min) / bin_size), 48) + 1;
bin(si == si_max) = 50;
p1 = zeros(1, 50);
for k = 1:50
    p1(k) = sum(si(bin == k)) / S;
end

% Compute weighted entropy
for I = 1:length(p1)
    pp(I) = (p1(I)*log(p1(I)));
end
pp(isnan(pp)) = 0;
w_ent = -1*(sum(pp));

end