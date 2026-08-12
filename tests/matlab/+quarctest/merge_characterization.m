function [S, per] = merge_characterization(S, per, S2, per2)
%MERGE_CHARACTERIZATION Fold a restricted pass into a full characterization.
%
%   [S, per] = quarctest.merge_characterization(S, per, S2, per2)
%
%   S, per    tables from a full run
%   S2, per2  tables from a later pass made with Metrics=..., adding columns
%             the full run did not carry
%
%   Returns the union, with rows ordered by system and then by the battery's
%   own metric order, so the report reads the same as if one run had produced
%   everything.
%
%   WHEN THIS IS LEGITIMATE. Only when the second pass used the same seed,
%   spread, R and N as the first. The ensemble is then the identical set of
%   initial conditions, the series generated from them are identical, and the
%   embedding is recomputed identically, so a metric measured in the second
%   pass belongs to exactly the realizations the first pass summarised. This
%   was verified before use rather than assumed: rerunning lyap_ros under a
%   restricted pass reproduced the full run's mean to 1.1e-16 on the logistic
%   map and 3.3e-16 on Lorenz, both at n = 1000.
%
%   The check is enforced here rather than left to the caller, because a merge
%   across differing protocols would silently produce a table whose rows came
%   from different ensembles while looking entirely ordinary.
%
%   Rows already present in S are kept; a metric appearing in both is taken
%   from S, so a merge cannot quietly rewrite an existing result.
%
%   See also QUARCTEST.CHARACTERIZE, QUARCTEST.WRITE_CHARACTERIZATION.

% Copyright (c) 2021-2026 Quantitative Analysis Research Core,
% Center for Human Movement Variability, University of Nebraska at Omaha.
% MIT licence. See LICENSE.txt.

arguments
    S   table
    per table
    S2  table
    per2 table
end

localAssertSameProtocol(S, S2);

% Keep only genuinely new (system, metric) pairs from the second pass, and drop
% its "not requested in this pass" placeholder rows.
S2 = S2(S2.role ~= "excluded" | S2.excludedBecause ~= "not requested in this pass", :);
haveKey = string(S.system) + "|" + string(S.metric);
newKey  = string(S2.system) + "|" + string(S2.metric);
S2 = S2(~ismember(newKey, haveKey), :);

S = [S; S2];

if ~isempty(per) && ~isempty(per2)
    haveKeyP = string(per.system) + "|" + string(per.metric);
    newKeyP  = string(per2.system) + "|" + string(per2.metric);
    per2 = per2(~ismember(newKeyP, unique(haveKeyP)), :);
    per = [per; per2];
end

% ---- restore the battery's ordering: catalogue order, then metric order

c = quarctest.sprott_catalog();
sysOrder = string([c.name]);
metOrder = string([quarctest.metric_policy(c(1)).id]);

si = localRank(string(S.system), sysOrder);
mi = localRank(string(S.metric), metOrder);
[~, ord] = sortrows([si mi]);
S = S(ord, :);
end

% ---------------------------------------------------------------- internals

function localAssertSameProtocol(S, S2)
f = {'N', 'R', 'seed', 'spread'};
for k = 1:numel(f)
    a = unique(S.(f{k})); b = unique(S2.(f{k}));
    if ~isequal(a, b)
        error('quarctest:merge:protocolMismatch', ...
            ['%s differs between the two passes (%s vs %s). The ensembles are ' ...
             'not the same realizations, so their rows cannot be merged.'], ...
            f{k}, mat2str(a), mat2str(b));
    end
end
end

function r = localRank(v, order)
r = zeros(numel(v), 1);
for i = 1:numel(order)
    r(v == order(i)) = i;
end
r(r == 0) = numel(order) + 1;
end
