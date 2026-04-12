function mu = computeReactionRates(obj, temperature)
    arguments
        obj  Model
        temperature double
    end
    mu = arrayfun(@(r) r.computeRate(temperature), obj.Reactions);
end