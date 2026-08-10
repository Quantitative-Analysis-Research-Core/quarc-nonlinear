function M = metric_policy(sys)
%METRIC_POLICY Which metrics are defensible on a given system, and why.
%
%   M = quarctest.metric_policy(sys)
%
%   Returns a struct array, one element per metric in the battery, with
%
%     id         column name used in the report
%     label      human-readable name
%     units      units of the value, or "" when dimensionless
%     applies    logical, whether this metric is computed for this system
%     reason     why it is excluded, when applies is false
%     reference  "lambda", "d2", or "" -- the catalogue field holding a
%                published value for this metric, if any
%
%   THE POINT OF THIS FUNCTION IS THE EXCLUSIONS. Every metric in the library
%   will return a number for any series handed to it. Whether that number
%   means anything depends on the system, and a report that quietly computes
%   all of them everywhere would put undefined quantities in a table next to
%   well-defined ones with no way to tell them apart. The exclusions are
%   recorded rather than applied silently, so a reader can disagree with a
%   judgement they can see.
%
%   ONLY TWO METRICS HAVE A PUBLISHED REFERENCE. Sprott's Appendix A gives the
%   largest Lyapunov exponent and the correlation dimension, the latter with a
%   stated uncertainty. Nothing in the appendix constrains a recurrence
%   determinism, a sample entropy, an AMI delay or an FNN dimension: those are
%   properties of an estimator applied to a series, not of the system, and no
%   reference value exists to compare them against. They are characterized --
%   measured across the ensemble and reported with their spread -- and they
%   are never validated. The report marks which is which on every row so the
%   distinction survives into whatever reads it.
%
%   See also QUARCTEST.METRIC_BATTERY, QUARCTEST.CHARACTERIZE.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    sys (1,1) struct
end

isMap = sys.kind == "map";

% Coarse-graining and scaling exponents both make a claim about behaviour
% ACROSS TIME SCALES. A map has no time scale: one iterate is the only unit
% there is, successive iterates are not samples of an underlying continuous
% trajectory, and averaging adjacent ones -- which is what multiscale
% coarse-graining does -- averages over an interval with no meaning. DFA's
% alpha is the same claim in a different form. Both are computed for flows,
% where the observable is a sampled continuous signal, and excluded for maps.
mapScaleReason = "coarse-graining and scaling exponents assume samples of " + ...
    "an underlying continuous-time signal; a map's iterates are not that";

M = struct('id', {}, 'label', {}, 'units', {}, 'applies', {}, ...
           'reason', {}, 'reference', {});

% AMI's delay and FNN's dimension are reported for every system, but on maps
% they do not drive the rest of the battery -- see quarctest.embed_policy.
% embed_delay and embed_dim record what was actually used, so a reader can
% always see both the estimate and the choice.
M = add(M, "ami_delay",  "AMI first-minimum delay", "samples", true,  "", "");
M = add(M, "fnn_dim",    "FNN embedding dimension", "",        true,  "", "");
M = add(M, "embed_delay","Delay used by the battery",     "samples", true, "", "");
M = add(M, "embed_dim",  "Dimension used by the battery", "",        true, "", "");

M = add(M, "lyap_wolf",  "Largest Lyapunov exponent (Wolf)", ...
         "nats per unit time", true, "", "lambda");
M = add(M, "lyap_ros",   "Largest Lyapunov exponent (Rosenstein)", ...
         "nats per unit time", true, "", "lambda");

% D2 is defined for the invariant measure of an attractor. Sprott publishes a
% value for the conservative systems too, so they are included and the
% category is carried into the report rather than the row being dropped.
M = add(M, "corr_dim",   "Correlation dimension", "", true, "", "d2");

M = add(M, "rqa_radius", "RQA radius at target recurrence", "", true, "", "");
M = add(M, "rqa_det",    "RQA determinism",   "percent", true, "", "");
M = add(M, "rqa_lam",    "RQA laminarity",    "percent", true, "", "");
M = add(M, "rqa_meanL",  "RQA mean diagonal line length", "samples", true, "", "");
M = add(M, "rqa_maxL",   "RQA longest diagonal line",     "samples", true, "", "");
M = add(M, "rqa_entL",   "RQA diagonal line entropy",     "bits", true, "", "");
M = add(M, "rqa_entV",   "RQA vertical line entropy",     "bits", true, "", "");
M = add(M, "rqa_entW",   "RQA weighted recurrence entropy", "bits", true, "", "");

M = add(M, "ent_samp",   "Sample entropy",      "", true, "", "");
M = add(M, "ent_ap",     "Approximate entropy", "", true, "", "");
M = add(M, "ent_permu",  "Permutation entropy", "bits", true, "", "");

M = add(M, "ent_ms_s1",  "Refined composite multiscale entropy, scale 1", "", ...
         ~isMap, ternary(isMap, mapScaleReason, ""), "");
M = add(M, "ent_ms_ci",  "Multiscale complexity index (sum over scales)", "", ...
         ~isMap, ternary(isMap, mapScaleReason, ""), "");
M = add(M, "dfa_alpha",  "DFA scaling exponent", "", ...
         ~isMap, ternary(isMap, mapScaleReason, ""), "");
end

% ---------------------------------------------------------------- internals

function M = add(M, id, label, units, applies, reason, reference)
M(end+1) = struct('id', string(id), 'label', string(label), ...
                  'units', string(units), 'applies', logical(applies), ...
                  'reason', string(reason), 'reference', string(reference)); %#ok<AGROW>
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
