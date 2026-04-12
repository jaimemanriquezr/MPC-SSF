function [K, denIdx, numIdx] = lookupQuotients(obj, components)
    names = [components.Name];
    dictHalfSaturationConstants = {obj.HalfSaturationConstants};
    keyComponents = cellfun(@(d) keys(d), dictHalfSaturationConstants, "UniformOutput", false);
    containsQuotient = cellfun(@(s) any(contains(s, "/")), keyComponents); 
    keyComponentsNum = cellfun(@(s, i) s + "/" + i, ...
                            keyComponents, num2cell(1:length(keyComponents)), "UniformOutput", false);
    keyQuotientsSplit = cellfun(@(s) split(s, "/", 2), ...
                            keyComponentsNum(containsQuotient), ...
                            "UniformOutput", false);

    keyQuotientsNames = cat(1, keyQuotientsSplit{:});
    [~, numIdx] = ismember(keyQuotientsNames(:, 1), names);
    [~, denIdx] = ismember(keyQuotientsNames(:, 2), names);
    rxId = double(keyQuotientsNames(:, 3));

    valueQuotients = cell2mat(cellfun(@(d, key) lookup(d, key), ...
                            dictHalfSaturationConstants(containsQuotient), ...
                            keyComponents(containsQuotient), ...
                            "UniformOutput", false));
                        
    K = nan(length(rxId), length(obj));
    for i = 1:length(rxId)
        K(i, rxId(i)) = valueQuotients(i);
    end
end