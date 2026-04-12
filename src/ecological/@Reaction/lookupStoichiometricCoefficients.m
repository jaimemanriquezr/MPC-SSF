function sigma = lookupStoichiometricCoefficients(obj, components)
    arguments
        obj   Reaction
        components  Component
    end
    names = [components.Name];

    numReactions = length(obj);
    numComponents = length(names);
    
    sigma = zeros(numComponents, numReactions);
    for i = 1:length(obj)
        sigma(:, i) = lookup(obj(i).StoichiometricCoefficients, names, FallbackValue=0);
    end
end