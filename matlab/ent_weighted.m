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
count = 1;
S = sum(si);
for s = si_min:bin_size:si_max
    P = sum(si(si>= s&si<(s+bin_size)));
    p1(count) = P / S;
    count = count+1;
end

% Compute weighted entropy
for I = 1:length(p1)
    pp(I) = (p1(I)*log(p1(I)));
end
pp(isnan(pp)) = 0;
w_ent = -1*(sum(pp));

end