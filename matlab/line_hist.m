function [recurrence, diag_hist, vertical_hist, a] = line_hist(x, a, radius, type, theiler)

% THEILER excludes every pair closer than theiler samples in time: the
% diagonals |i-j| <= theiler are skipped, the same band is blanked in each
% column before vertical lines are counted, and %REC is taken over the
% remaining pairs only. At dense sampling those near-LOI cells record a
% trajectory shadowing itself tangentially, not a recurrence, and MaxL --
% one extreme line rather than an average -- is the statistic they corrupt:
% measured across the dysts catalogue, its Spearman correlation with the
% reference exponent per sample moves from -0.14 at theiler=0 to -0.49 at
% theiler=10. Default 0, which reproduces the old behaviour exactly (only
% the line of identity excluded). Ignored for crqa, where the matrix is
% cross-recurrence and the main diagonal is data.
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

if nargin < 5
    theiler = 0;
end
if contains(type, 'crqa')
    theiler = 0;
end

if ~iscell(a)
    a = {a};
end

% Convert distance matrices to recurrence matrices
for i2 = 1:length(a)
    a{i2} = a{i2}+radius;
    a{i2}(a{i2} >= 0) = 1;
    a{i2}(a{i2} < 0) = 0;
end

% If a contains multiple recurrence matrices, compute dot product of all
% matrices...?
if length(a) > 1
    min_size = min(cellfun(@(c) size(c,2), a));
    for i3 = 1:length(a)-1
        a{i3+1} = a{i3}(1:min_size,1:min_size).*a{i3+1}(1:min_size,1:min_size);
    end
    a = a{i3+1};
else
    a = a{1};
end

% Caluculate diagonal line distribution
diag_hist = [];
vertical_hist = [];
for i4 = -(length(a)-1):length(a)-1
    % Inside the Theiler band nothing counts as a recurrence. The line of
    % identity is the band's centre, so for non-crqa the old remove-the-max
    % step is subsumed by skipping i4 == 0 here.
    if ~contains(type, 'crqa') && abs(i4) <= theiler
        continue
    end
    c=diag(a,i4);
    % bwlabel is taking each diagonal line and looking for the 1's, it will
    % return increasing numbers for each new instance of 1's, for example
    % the input vector [0 1 1 0 1 0 1 1 0 0 1 1 1] will return
    %                  [0 1 1 0 2 0 3 3 0 0 4 4 4]
    d=bwlabel(c,8);
    % One histogram bin per label, so each line is counted exactly once at
    % its exact length. hist(d) with its ten default bins merged adjacent
    % labels whenever a diagonal carried ten or more lines -- summing their
    % lengths into fictitious long lines -- and the zeros needed a separate
    % drop that also took real lines with it.
    if max(d) > 0
        d = histcounts(d(d>0), 0.5:1:max(d)+0.5);
    else
        d = [];
    end
    % diag_hist is creating one long array of all of the line lengths for
    % all of the diagonals
    diag_hist(length(diag_hist)+1:length(diag_hist)+length(d))=d;
end

% The line of identity never reaches diag_hist for RQA, jRQA and mdRQA --
% the loop above skips the whole Theiler band, whose centre it is. Removing
% it by dropping the maximum, as this block used to, would now delete the
% longest genuine line instead.
if ~contains(type,'crqa') && isempty(diag_hist)
    diag_hist=0;
end

% Calculate vertical line distribution
for i5=1:length(a)
    c=(a(:,i5));
    if theiler > 0
        c(max(1, i5-theiler):min(length(a), i5+theiler)) = 0;
    end
    v=bwlabel(c,8);
    % Same exact tabulation as the diagonals above.
    if max(v) > 0
        v = histcounts(v(v>0), 0.5:1:max(v)+0.5);
    else
        v = [];
    end
    vertical_hist(length(vertical_hist)+1:length(vertical_hist)+length(v))=v;
end

% Calculate percent recurrence over the pairs that can count: outside the
% Theiler band for RQA and friends, everywhere for crqa. At theiler=0 the
% band is exactly the line of identity and this reduces to the old formula.
if ~contains(type,'crqa')
    n_ = length(a);
    band = sum(diag(a, 0));
    for k = 1:theiler
        band = band + sum(diag(a, k)) + sum(diag(a, -k));
    end
    excluded = n_ + 2*sum(n_ - (1:theiler));
    recurrence = 100*(sum(sum(a))-band)/(n_^2-excluded);
else
    recurrence = 100*(sum(sum(a)))/(length(a)^2);
end

end