function E = embed_policy(sys)
%EMBED_POLICY The embedding used to drive the battery, and why.
%
%   E = quarctest.embed_policy(sys)
%
%   Fields
%     delayRule  "fixed" or "ami"
%     delayValue the delay, when delayRule is "fixed"
%     dimRule    "fixed" or "fnn"
%     dimValue   the dimension, when dimRule is "fixed"
%     stateDim   the system's state dimension, from the catalogue
%     reason     one sentence on why this system gets this rule
%
%   THE EMBEDDING IS A PROTOCOL, NOT A MEASUREMENT. Every metric downstream of
%   the reconstruction -- correlation dimension, both Lyapunov exponents, the
%   whole RQA block -- depends on the delay and dimension it is given, often
%   more strongly than on the estimator being characterized. So the embedding
%   is stated here, per system, rather than being whatever the library's
%   estimators happened to return on that realization.
%
%   AMI's delay and FNN's dimension are still computed and reported as metrics
%   in their own right. What changes is that on maps they do not drive the
%   rest of the battery. Keeping the two roles separate means the report can
%   show that FNN returns 1 on the logistic map and 5 on Ikeda without those
%   answers silently corrupting every other column.
%
%   MAPS: DELAY 1, BY CONSTRUCTION. A map's successive iterates are its
%   natural coordinates. There is no sampling interval to choose, and the
%   first minimum of AMI on a map is not an embedding delay -- it is a number
%   produced by treating a recursion as though it were a sampled continuous
%   signal. Measured on the logistic map, taking AMI's answer of 6 with FNN's
%   dimension yields a correlation dimension of 3.09 against a published D2 of
%   exactly 1.0; at delay 1 the same estimator returns 0.93.
%
%   MAPS: DIMENSION ONE ABOVE THE STATE. The catalogue knows each map's state
%   dimension, so the reconstruction does not have to guess it. One coordinate
%   above the state suffices to unfold a map observed through a single
%   variable and its history. This matches the convention already used in
%   quarctest.lye_benchmark (2 for one-dimensional maps, 3 for two-dimensional
%   ones). Measured against Sprott's D2: logistic 0.931 vs 1.000, cusp 0.999
%   vs 1.000, Henon 1.219 vs 1.220, Ikeda 1.691 vs 1.690.
%
%   FLOWS: THE LIBRARY'S OWN ESTIMATORS. For a flow sampled through one
%   observable, the delay and dimension are genuinely unknown and choosing
%   them is part of the analysis -- which is the situation a user of this
%   library is actually in with a recording. So flows use ami() and fnn(),
%   and the scatter those choices induce across the ensemble is part of what
%   the report is measuring. The state dimension is used only as a floor, to
%   stop a failed FNN search from collapsing the reconstruction.
%
%   See also QUARCTEST.METRIC_POLICY, QUARCTEST.METRIC_BATTERY.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    sys (1,1) struct
end

d = numel(sys.x0);

if sys.kind == "map"
    E = struct('delayRule', "fixed", 'delayValue', 1, ...
               'dimRule', "fixed", 'dimValue', d + 1, ...
               'stateDim', d, ...
               'reason', "map: iterates are the natural coordinates, so the " + ...
                         "delay is 1 by construction and the state dimension " + ...
                         "is known rather than estimated");
else
    E = struct('delayRule', "ami", 'delayValue', NaN, ...
               'dimRule', "fnn", 'dimValue', NaN, ...
               'stateDim', d, ...
               'reason', "flow: delay and dimension are unknown for a scalar " + ...
                         "observable and are estimated the way a user would " + ...
                         "estimate them, so their scatter is part of the result");
end
end
