function epsilon = computePorosity(obj, z)
    arguments
        obj  SandFilter
        z double
    end
    epsilon0 = obj.SandPorosity;
    delta = obj.SandRoughness;
    epsilon = bound((epsilon0 - 1)/delta*z + epsilon0, [epsilon0 1]);
end

function m = bound(x, bounds)
    m = max(bounds(1), min(x, bounds(2)));
end