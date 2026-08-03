function [RP, RESULTS]=jrqa(data,tau,dim,param,threshold,options)
%JRQA Joint recurrence quantification analysis of several time series.
%   [RP,RESULTS] = JRQA(DATA,TAU,DIM) reconstructs a phase space from each
%   column of DATA using that channel's own delay TAU(i) and embedding
%   dimension DIM(i), then returns the joint recurrence plot RP and a
%   struct RESULTS of recurrence variables.
%
%   DATA may instead already be a cell array of per-channel phase spaces,
%   {Y1, Y2, ...}, when options.PhaseSpace is true: see below. Joint
%   recurrence is inherently per-channel, so unlike RQA, CRQA and MDRQA a
%   single N-by-D matrix cannot unambiguously stand in for a phase space:
%   the channel boundaries would be lost.
%
%   ___ = JRQA(DATA,TAU,DIM,PARAM,THRESHOLD) selects what THRESHOLD means:
%      PARAM = "rad"  THRESHOLD is the recurrence radius directly.
%      PARAM = "rec"  (default) THRESHOLD is a target percent recurrence;
%                      the radius is searched for iteratively.
%
%   ___ = JRQA(...,Name=Value) sets additional options:
%      Zscore       (1,1) 0 or 1. Z-score each channel before use. Default 1.
%      Norm         "euc", "max", "min" or "none". Distance-matrix
%                   normalization. Default "none".
%      Dmin, Vmin   Minimum diagonal/vertical line length counted as
%                   deterministic/laminar. Default 2.
%      Plot         (1,1) 0 or 1. Show the recurrence plot. Default 0.
%      Iter         Bisection iterations used to search for the radius
%                   under PARAM="rec". Default 20.
%      PhaseSpace   (1,1) logical. If true, DATA is a cell array with one
%                   already-reconstructed phase space per channel and is
%                   used verbatim: TAU and DIM are recorded in RESULTS but
%                   no reconstruction is performed. Default false.
%
%   Notes
%   A supplied phase space carries no record of the delay or dimension used
%   to build it, so Norm="euc" -- whose scaling depends on DIM matching each
%   phase space's actual width -- cannot verify that assumption and JRQA
%   warns rather than silently trusting a possibly stale DIM.
%
%   Examples
%      % Reconstruct internally
%      [rp, r] = jrqa([x y], tau, dim);
%
%      % Reconstruct once per channel, share it, skip re-embedding
%      Y1 = psr(x, tau(1), dim(1));
%      Y2 = psr(y, tau(2), dim(2));
%      [rp, r] = jrqa({Y1, Y2}, tau, dim, PhaseSpace=true);
%
%   See also PSR, RQA, CRQA, MDRQA.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    data {mustBeA(data, ["double","cell"])}
    tau (1,2) {mustBeInteger, mustBePositive} = [1 1]
    dim (1,2) {mustBeInteger, mustBePositive} = [1 1]
    param (1,1) string {mustBeMember(param,["rad", "rec"])} = "rec"
    threshold (1,1) double {mustBePositive} = 2.5
    options.Zscore (1,1) {mustBeMember(options.Zscore,[0,1])} = 1
    options.Norm (1,1) {mustBeMember(options.Norm,["euc", "max", "min", "none"])} = "none"
    options.Dmin (1,1) {mustBeInteger, mustBePositive} = 2
    options.Vmin (1,1) {mustBeInteger, mustBePositive} = 2
    options.Plot (1,1) {mustBeMember(options.Plot,[0,1])} = 0
    options.Orient (1,1) {mustBeMember(options.Orient,["col", "row"])} = "col"
    options.Iter (1,1) {mustBeInteger, mustBePositive} = 20
    options.PhaseSpace (1,1) logical = false
end

% options is not visible yet when the arguments block validates data, so
% the type/shape checks that depend on it live here instead.
if options.PhaseSpace
    if ~iscell(data)
        error('jrqa:phaseSpaceMustBeCell', ...
            ['PhaseSpace is true, so DATA must be a cell array with one ' ...
             'already-reconstructed phase space per channel, e.g. ' ...
             '{Y1, Y2}, rather than a raw multi-column series.']);
    elseif numel(data) < 2
        error('jrqa:tooFewChannels', ...
            'DATA must contain at least two channels.');
    end
elseif ~isa(data, 'double')
    error('jrqa:rawDataMustBeDouble', ...
        'DATA must be a double matrix unless PhaseSpace is true.');
else
    mustbeAtLeastTwoColumns(data)
end

%% Begin code

%% Change variable names for readability
dmin = options.Dmin;
vmin = options.Vmin;

%% Standardize data if zscore is true
% If zscore is selected then zscore the data (each channel independently
% when a cell array of phase spaces was supplied)
if options.Zscore
    if options.PhaseSpace
        data = cellfun(@zscore, data, 'UniformOutput', false);
    else
        data = zscore(data);
    end
end

if options.PhaseSpace
    % A supplied phase space is used verbatim; TAU and DIM cannot be
    % confirmed to match it.
    data2 = data;
    DIM = numel(data2);
    if options.Norm == "euc"
        warning('jrqa:eucNormAssumesDim', ...
            ['PhaseSpace is true, so DATA was not rebuilt from TAU and DIM ' ...
             'and DIM cannot be confirmed to match its width. Norm="euc" ' ...
             'scales by DIM regardless; pass DIM to match each phase ' ...
             'space''s width or use Norm="none", "min" or "max".']);
    end
else
    % Get number of time series
    DIM = size(data, 2);

    % Embed the data onto phase space
    for i = 1:DIM
        if dim(i) > 1
            data2{i} = psr(data(:,i), tau(i), dim(i));
        end
    end
end

% Calculate distance matrix based on the type of RQA
for i = 1:DIM
    a{i}=pdist2(data2{i},data2{i});
    a{i}=abs(a{i})*-1;
end

% Normalize distance matrix
if contains(options.Norm, 'euc')
    for i = 1:length(a)
        b = mean(a{i}(a{i}<0));
        b = -sqrt(abs(((b^2)+2*(DIM*dim))));
        a{i} = a{i}/abs(b);
    end
elseif contains(options.Norm, 'min')
    for i = 1:length(a)
        b = max(a{i}(a{i}<0));
        a{i} = a{i}/abs(b);
    end
elseif contains(options.Norm, 'max')
    for i = 1:length(a)
        b = min(a{i}(a{i}<0));
        a{i} = a{i}/abs(b);
    end
end

% Compute weighted recurrence plot
wrp = a{1};
% for i = 1:size(a,2)-1
%     wrp{i+1} = wrp{i}.*wrp{i+1};
% end
% if i
%     wrp = -(abs(wrp{i+1})).^(1/(i+1));
% end
% if iscell(wrp)
%     wrp = wrp{1};
% end

% Calculate recurrence plot
switch param
    case 'rad'
        % THRESHOLD is the radius itself in this branch.
        radius = threshold;
        [recurrence, diag_hist, vertical_hist,A] = line_hist(data,a,threshold,'jrqa');
    case 'rec'
        radius_start = 0.01;
        radius_end = 0.5;
        [recurrence, diag_hist, vertical_hist, radius, A] = set_radius(data,a,radius_start,radius_end,threshold,'jrqa',options.Iter);
end

%% Calculate RQA variabes
RESULTS.DIM = 1;
RESULTS.EMB = dim;
RESULTS.DEL = tau;
RESULTS.RADIUS = radius;
RESULTS.NORM = options.Norm;
RESULTS.ZSCORE = options.Zscore;
RESULTS.Size=length(A);
RESULTS.REC = recurrence;
if RESULTS.REC > 0
    RESULTS.DET=100*sum(diag_hist(diag_hist>=dmin))/sum(diag_hist);
    RESULTS.MeanL=mean(diag_hist(diag_hist>=dmin));
    RESULTS.MaxL=max(diag_hist(diag_hist>=dmin));
    [count,bin]=hist(diag_hist(diag_hist>=dmin),min(diag_hist(diag_hist>=dmin)):max(diag_hist(diag_hist>=dmin)));
    total=sum(count);
    p=count./total;
    del=find(count==0); p(del)=[];
    RESULTS.EntrL=-sum(p.*log2(p));
    RESULTS.LAM=100*sum(vertical_hist(vertical_hist>=vmin))/sum(vertical_hist);
    RESULTS.MeanV=mean(vertical_hist(vertical_hist>=vmin));
    RESULTS.MaxV=max(vertical_hist(vertical_hist>=vmin));
    [count,bin]=hist(vertical_hist(vertical_hist>=vmin),min(vertical_hist(vertical_hist>=vmin)):max(vertical_hist(vertical_hist>=vmin)));
    total=sum(count);
    p=count./total;
    del=find(count==0); p(del)=[];
    RESULTS.EntrV=-sum(p.*log2(p));
    RESULTS.EntrW=NaN;
else
    RESULTS.DET=NaN;
    RESULTS.MeanL=NaN;
    RESULTS.MaxL=NaN;
    RESULTS.EntrL=NaN;
    RESULTS.LAM=NaN;
    RESULTS.MeanV=NaN;
    RESULTS.MaxV=NaN;
    RESULTS.EntrV=NaN;
    RESULTS.EntrW=NaN;
end
RP=imrotate(1-A,90);

%% Plot
if options.Plot
    rqa_plot(data, RP, RESULTS, tau, dim, DIM, options.Zscore, options.Norm, radius, wrp, 'jrqa');
end

end

% Custom validation function
function mustbeAtLeastTwoColumns(data)
% Test for size
if size(data,2) < 2
    error('Data must be at least two column vectors.')
end
end