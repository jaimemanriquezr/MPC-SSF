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
    % theta^(T/T_nom - 1). The ratio form appears in ecomodel.tex as a typo
    % (would make the response ~293x too weak); the legacy run_biofilm.m:91
    % uses theta^(T_K - 293) -- the SAME difference form as here (an earlier
    % version of this comment claimed it was sign-flipped; verified wrong
    % 2026-08-20).
    switch obj.TemperatureResponse
        case "exponential"
            mu = mu20 .* (theta).^(temperature - nominalTemperature);
        case "cardinal"
            % CTMI (Rosso et al. 1993), evaluated in degC, rescaled so that
            % mu(nominalTemperature) = NominalRate exactly (pre-registered
            % anchor: Table 2 comparability, 19-20 degC behavior preserved).
            tC = temperature - 273;
            tRef = nominalTemperature - 273;
            ct = obj.CardinalTemperatures;
            mu = mu20 * ctmi(tC, ct) / ctmi(tRef, ct);
        otherwise
            error("Invalid TemperatureResponse.")
    end
end

function phi = ctmi(t, ct)
    tmin = ct(1); topt = ct(2); tmax = ct(3);
    if t <= tmin || t >= tmax
        phi = 0.0;
        return
    end
    num = (t - tmax) * (t - tmin)^2;
    den = (topt - tmin) * ((topt - tmin)*(t - topt) - (topt - tmax)*(topt + tmin - 2*t));
    phi = num / den;
end