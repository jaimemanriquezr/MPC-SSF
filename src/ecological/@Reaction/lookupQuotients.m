function [K, denIdx, numIdx] = lookupQuotients(obj, components)
    names = [components.Name];
    dictHalfSaturationConstants = {obj.HalfSaturationConstants};
    keyComponents = cellfun(@(d) keys(d), dictHalfSaturationConstants, "UniformOutput", false);
    containsQuotient = cellfun(@(s) any(contains(s, "/")), keyComponents); 
    % A reaction may mix plain and quotient half-saturations (e.g. O2 alongside
    % PG/PHO): keep only the quotient keys before tagging with the reaction index.
    keyComponentsQ = cellfun(@(s) s(contains(s, "/")), keyComponents, "UniformOutput", false);
    keyComponentsNum = cellfun(@(s, i) s + "/" + i, ...
                            keyComponentsQ, num2cell(1:length(keyComponentsQ)), "UniformOutput", false);
    keyQuotientsSplit = cellfun(@(s) split(s, "/", 2), ...
                            keyComponentsNum(containsQuotient), ...
                            "UniformOutput", false);

    % No quotient half-saturations anywhere: return empty (0 x nRx) so the
    % [halfSaturationK; quotientK] concatenation in simulate degrades cleanly.
    if isempty(keyQuotientsSplit)
        K = nan(0, length(obj));  denIdx = [];  numIdx = [];
        return
    end
    keyQuotientsNames = cat(1, keyQuotientsSplit{:});
    [~, numIdx] = ismember(keyQuotientsNames(:, 1), names);
    [~, denIdx] = ismember(keyQuotientsNames(:, 2), names);
    rxId = double(keyQuotientsNames(:, 3));

    valueQuotients = cell2mat(cellfun(@(d, key) lookup(d, key), ...
                            dictHalfSaturationConstants(containsQuotient), ...
                            keyComponentsQ(containsQuotient), ...
                            "UniformOutput", false));
                        
    K = nan(length(rxId), length(obj));
    for i = 1:length(rxId)
        K(i, rxId(i)) = valueQuotients(i);
    end
end