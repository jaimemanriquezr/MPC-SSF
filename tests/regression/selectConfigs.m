function sel = selectConfigs(cfgs, names)
% SELECTCONFIGS  Indices of the configs matching `names` (exact or wildcard).
%
%   sel = selectConfigs(regressionConfigs(), ["theta_*", "film_transfer"])
%
% Errors on a pattern that matches nothing, rather than silently running an
% empty suite -- a green "0 passed, 0 failed" is worse than a red one.

all_ = string({cfgs.name});
if isscalar(names) && names == "*"
    sel = 1:numel(cfgs);
    return
end

hit = false(1, numel(cfgs));
for p = names
    m = ~cellfun(@isempty, regexp(all_, "^" + regexptranslate("wildcard", p) + "$", "once"));
    if ~any(m)
        error("selectConfigs:noMatch", ...
            "'%s' matches no configuration. Available: %s", p, strjoin(all_, ", "));
    end
    hit = hit | m;
end
sel = find(hit);
end