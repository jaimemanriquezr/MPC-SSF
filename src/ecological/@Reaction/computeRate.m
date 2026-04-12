function mu = computeRate(obj, temperature, options)
    arguments
        obj   Reaction
        temperature double
        options.TemperatureScale = "Celsius";
        options.NominalTemperature = 20;
    end
    mu20 = obj.NominalRate;
    theta = obj.TemperatureCorrectionFactor;

    switch options.TemperatureScale
        case "Celsius"
            temperature = 273 + temperature;
            nominalTemperature = 273 + options.NominalTemperature;
        case "Kelvin"
            nominalTemperature = options.NominalTemperature;
        otherwise
            error("Invalid temperature scale.") 
    end
    mu = mu20 .* (theta).^(temperature/nominalTemperature - 1);
end