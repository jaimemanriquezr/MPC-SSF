function p = lookupOrder(obj, components)
    arguments
        obj mpcssf.Reaction
        components mpcssf.Component
    end
    names = [components.Name];

    numReactions = length(obj);
    numComponents = length(names);
    
    p = zeros(numComponents, numReactions);
    for i = 1:length(obj)
        p(:, i) = lookup(obj(i).Order, names, FallbackValue=0);
    end
end