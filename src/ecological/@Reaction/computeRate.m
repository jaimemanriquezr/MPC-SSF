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
    % theta (~1.05-1.08) is calibrated to the temperature DIFFERENCE form
    % theta^(T - T_nom) (a K/degC difference), not the dimensionless ratio
    % theta^(T/T_nom - 1). The ratio makes the response ~293x too weak unless
    % theta is recalibrated to theta^293; kept the calibrated theta and the
    % difference form. (Authoritative slow-sand run_biofilm uses theta^(293-T_K),
    % sign-flipped; all three coincide only at the 20 degC reference.)
    mu = mu20 .* (theta).^(temperature - nominalTemperature);
end