function [RP, RESULTS]=crqa(x,delay,dim,param,threshold,options)
%CRQA Cross recurrence quantification analysis of two time series.
%   [RP,RESULTS] = CRQA(DATA,TAU,DIM) reconstructs a phase space from each
%   column of the two-column matrix DATA using delay TAU and embedding
%   dimension DIM, then returns the cross recurrence plot RP and a struct
%   RESULTS of recurrence variables.
%
%   DATA may instead already be a phase space (an N-by-2*DIM matrix, with
%   the two series' lagged coordinates interleaved as PSR produces them)
%   when options.phasespace is true: see below.
%
%   ___ = CRQA(DATA,TAU,DIM,PARAM,THRESHOLD) selects what THRESHOLD means:
%      PARAM = "rad"  THRESHOLD is the recurrence radius directly.
%      PARAM = "rec"  (default) THRESHOLD is a target percent recurrence;
%                      the radius is searched for iteratively.
%
%   ___ = CRQA(...,Name=Value) sets additional options:
%      Zscore       (1,1) 0 or 1. Z-score DATA before use. Default 1.
%      Norm         "euc", "max", "min" or "none". Distance-matrix
%                   normalization. Default "none".
%      dmin, Vmin   Minimum diagonal/vertical line length counted as
%                   deterministic/laminar. Default 2.
%      Plot         (1,1) 0 or 1. Show the recurrence plot. Default 0.
%      Iter         Bisection iterations used to search for the radius
%                   under PARAM="rec". Default 20.
%      PhaseSpace   (1,1) logical. If true, DATA is already a
%                   reconstructed phase space and is used verbatim: TAU
%                   and DIM are recorded in RESULTS but no reconstruction
%                   is performed. Default false.
%
%   Notes
%   A supplied phase space carries no record of the delay or dimension used
%   to build it, so norm="euc" -- whose scaling depends on DIM matching the
%   phase space's actual width -- cannot verify that assumption and CRQA
%   warns rather than silently trusting a possibly stale DIM.
%
%   Examples
%      % Reconstruct internally
%      [rp, r] = crqa([x y], delay, dim);
%
%      % Reconstruct once, share it, skip re-embedding
%      Y = psr([x y], delay, dim);
%      [rp, r] = crqa(Y, delay, dim, phasespace=true);
%
%   See also PSR, RQA, JRQA, MDRQA.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    x double {mustBeNonempty}
    delay (1,1) {mustBeInteger, mustBePositive} = 1
    dim (1,1) {mustBeInteger, mustBePositive} = 1
    param (1,1) string {mustBeMember(param,["rad", "rec"])} = "rec"
    threshold (1,1) double {mustBePositive} = 2.5
    options.zscore (1,1) {mustBeMember(options.zscore,[0,1])} = 1
    options.norm (1,1) {mustBeMember(options.norm,["euc", "max", "min", "none"])} = "none"
    options.dmin (1,1) {mustBeInteger, mustBePositive} = 2
    options.vmin (1,1) {mustBeInteger, mustBePositive} = 2
    options.plot (1,1) {mustBeMember(options.plot,[0,1])} = 0
    options.orient (1,1) {mustBeMember(options.orient,["col", "row"])} = "col"
    options.iter (1,1) {mustBeInteger, mustBePositive} = 20
    options.phasespace (1,1) logical = false
end

% mustBeTwoColumns only applies to a raw pair of series: a supplied phase
% space is legitimately wider. options is not visible yet when the
% arguments block validates x, so the shape is checked here instead.
if ~options.phasespace
    mustBeTwoColumns(x)
end

%% Begin code

%% Change variable names for readability
dmin = options.dmin;
vmin = options.vmin;

%% Standardize x if zscore is true
% If zscore is selected then zscore the x
if options.zscore
    x = zscore(x);
end

% A supplied phase space is used verbatim; otherwise embed the x onto
% phase space as before.
if options.phasespace
    if options.norm == "euc"
        warning('crqa:eucNormAssumesDim', ...
            ['PhaseSpace is true, so DATA was not rebuilt from TAU and DIM ' ...
             'and DIM cannot be confirmed to match its width. Norm="euc" ' ...
             'scales by DIM regardless; pass DIM equal to size(data,2)/2 or ' ...
             'use Norm="none", "min" or "max".']);
    end
elseif dim > 1
    x = psr(x, delay, dim);
end

% Calculate distance matrix based on the type of RQA
a = pdist2(x(:,1:2:end),x(:,2:2:end));
a = abs(a)*-1;

% Normalize distance matrix
if contains(options.norm, 'euc')
    b = mean(a(a<0));
    b = -sqrt(abs(((b^2)+2*(2*dim))));
    a = a/abs(b);
elseif contains(options.norm, 'min')
    b = max(a(a<0));
    a = a/abs(b);
elseif contains(options.norm, 'max')
    b = min(a(a<0));
    a = a/abs(b);
end

% % Compute weighted recurrence plot (doesn't seem like it's doing anything
% % for univariate RQA, cross RQA)
% wrp = a;
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
        [recurrence, diag_hist, vertical_hist,A] = line_hist(x,a,threshold,'crqa');
    case 'rec'
        radius_start = 0.01;
        radius_end = 0.5;
        [recurrence, diag_hist, vertical_hist, radius, A] = set_radius(x,a,radius_start,radius_end,threshold,'crqa',options.iter);
end

%% Calculate RQA variabes
RESULTS.DIM = 1;
RESULTS.EMB = dim;
RESULTS.DEL = delay;
RESULTS.RADIUS = radius;
RESULTS.NORM = options.norm;
RESULTS.ZSCORE = options.zscore;
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
    RESULTS.EntrW=ent_weighted(a);
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
if options.plot
    rqa_plot(x, RP, RESULTS, delay, dim, 2, options.zscore, options.norm, radius, a, 'crqa');
end

end

% Custom validation function
function mustBeTwoColumns(data)
% Test for size
if size(data,2) ~= 2
    error('Data must be two column vectors.')
end
end