function [recurrence, diag_hist, vertical_hist, a] = line_hist(x, a, radius, type)

% Check if a is cell
% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

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

% Remove the line of identity in RQA, jRQA, and mdRQA
if ~contains(type,'crqa')
    diag_hist=diag_hist(diag_hist<max(diag_hist));
    if isempty(diag_hist)
        diag_hist=0;
    end
end

% Calculate vertical line distribution
for i5=1:length(a)
    c=(a(:,i5));
    v=bwlabel(c,8);
    % Same exact tabulation as the diagonals above.
    if max(v) > 0
        v = histcounts(v(v>0), 0.5:1:max(v)+0.5);
    else
        v = [];
    end
    vertical_hist(length(vertical_hist)+1:length(vertical_hist)+length(v))=v;
end

% Calculate percent recurrence
if ~contains(type,'crqa')
    recurrence = 100*(sum(sum(a))-length(a))/(length(a)^2-length(a));
else
    recurrence = 100*(sum(sum(a)))/(length(a)^2);
end

end