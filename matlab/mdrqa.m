function [RP, RESULTS]=mdrqa(data,tau,dim,param,threshold,options)
%MDRQA Multidimensional recurrence quantification analysis.
%   [RP,RESULTS] = MDRQA(DATA,TAU,DIM) reconstructs a phase space from the
%   (possibly multi-column) matrix DATA using delay TAU and embedding
%   dimension DIM, applied uniformly across all columns, then returns the
%   recurrence plot RP and a struct RESULTS of recurrence variables. With
%   DIM=1 (the default) DATA is used as a phase space directly, with no
%   reconstruction: this is the usual way to feed MDRQA several already
%   co-registered channels (e.g. joint angles) with no added lag embedding.
%
%   DATA may also already be a phase space built with DIM>1 when
%   options.PhaseSpace is true: see below.
%
%   ___ = MDRQA(DATA,TAU,DIM,PARAM,THRESHOLD) selects what THRESHOLD means:
%      PARAM = "rad"  THRESHOLD is the recurrence radius directly.
%      PARAM = "rec"  (default) THRESHOLD is a target percent recurrence;
%                      the radius is searched for iteratively.
%
%   ___ = MDRQA(...,Name=Value) sets additional options:
%      Zscore       (1,1) 0 or 1. Z-score DATA before use. Default 1.
%      Norm         "euc", "max", "min" or "none". Distance-matrix
%                   normalization. Default "none".
%      Dmin, Vmin   Minimum diagonal/vertical line length counted as
%                   deterministic/laminar. Default 2.
%      Plot         (1,1) 0 or 1. Show the recurrence plot. Default 0.
%      Iter         Bisection iterations used to search for the radius
%                   under PARAM="rec". Default 20.
%      PhaseSpace   (1,1) logical. If true, DIM>1 is not used to trigger
%                   reconstruction: DATA is already a reconstructed phase
%                   space and is used verbatim. TAU and DIM are still
%                   recorded in RESULTS. Default false.
%
%   Notes
%   Because DATA can legitimately be multi-column either as raw multichannel
%   input or as an already-built phase space, PhaseSpace is what decides
%   which one it is here -- width alone is ambiguous. A supplied phase space
%   carries no record of the delay or dimension used to build it, so
%   Norm="euc" -- whose scaling depends on DIM matching the phase space's
%   actual width -- cannot verify that assumption and MDRQA warns rather
%   than silently trusting a possibly stale DIM.
%
%   Examples
%      % Several raw channels used directly, no lag embedding
%      [rp, r] = mdrqa([x y z]);
%
%      % Reconstruct once, share it, skip re-embedding
%      Y = psr([x y z], tau, dim);
%      [rp, r] = mdrqa(Y, tau, dim, PhaseSpace=true);
%
%   References
%      Wallot, S., Roepstorff, A. and Monster, D. (2016). Multidimensional
%      recurrence quantification analysis (MdRQA) for the analysis of
%      multidimensional time-series: A software implementation in MATLAB
%      and its application to group-level data in joint action. Frontiers
%      in Psychology, 7, 1835.
%
%   See also PSR, RQA, CRQA, JRQA.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    data double
    tau (1,1) {mustBeInteger, mustBePositive} = 1
    dim (1,1) {mustBeInteger, mustBePositive} = 1
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

%% Begin code

%% Change variable names for readability
dmin = options.Dmin;
vmin = options.Vmin;

%% Standardize data if zscore is true
% If zscore is selected then zscore the data
if options.Zscore
    data = zscore(data);
end

% A supplied phase space is used verbatim; otherwise embed the data onto
% phase space as before. DATA can legitimately be multi-column either way
% (several raw channels, or an already-built phase space), so PhaseSpace,
% not width, is what decides whether DIM>1 triggers reconstruction.
if options.PhaseSpace
    if options.Norm == "euc"
        warning('mdrqa:eucNormAssumesDim', ...
            ['PhaseSpace is true, so DATA was not rebuilt from TAU and DIM ' ...
             'and DIM cannot be confirmed to match its width. Norm="euc" ' ...
             'scales by DIM regardless; pass DIM equal to size(data,2) or ' ...
             'use Norm="none", "min" or "max".']);
    end
elseif dim > 1
    data = psr(data, dim, tau);
end

% Calculate distance matrix based on the type of RQA
a = pdist2(data,data);
a = abs(a)*-1;

% Normalize distance matrix
if contains(options.Norm, 'euc')
        b = mean(a(a<0));
        b = -sqrt(abs(((b^2)+2*(1*dim))));
        a = a/abs(b);
elseif contains(options.Norm, 'min')
        b = max(a(a<0));
        a = a/abs(b);
elseif contains(options.Norm, 'max')
        b = min(a(a<0));
        a = a/abs(b);
end

% % Compute weighted recurrence plot (doesn't seem like it's doing anything
% % for univariate RQA, cross RQA, multidimensional RQA)
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
        [recurrence, diag_hist, vertical_hist,A] = line_hist(data,a,threshold,'mdrqa');
    case 'rec'
        radius_start = 0.01;
        radius_end = 0.5;
        [recurrence, diag_hist, vertical_hist, radius, A] = set_radius(data,a,radius_start,radius_end,threshold,'mdrqa',options.Iter);
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
if options.Plot
    rqa_plot(data, RP, RESULTS, tau, dim, 1, options.Zscore, options.Norm, radius, a, 'mdrqa');
end

end

% Custom validation function
function mustBeSingleColumn(data)
    % Test for size
    if size(data,2) ~= 1
        error('Data must be single column vector.')
    end
end