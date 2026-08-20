% Cardinal temperature response (CTMI, Rosso 1993) unit checks - mirrors the
% Julia testset "cardinal temperature response".
W = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(W, "src")));
ct = [0.0, 25.0, 35.0];
rx = Reaction(NominalRate=3.0, TemperatureResponse="cardinal", ...
    CardinalTemperatures=ct);
assert(abs(rx.computeRate(20) - 3.0) < 1e-12, "anchor mu(20) failed");
assert(rx.computeRate(0) == 0.0 && rx.computeRate(-5) == 0.0 && ...
    rx.computeRate(35) == 0.0, "cardinal bounds failed");
assert(rx.computeRate(25) > 3.0, "Topt peak failed");
phi = @(t) ((t - 35).*t.^2) ./ (25*(25*(t - 25) - (-10)*(25 - 2*t)));
assert(abs(rx.computeRate(3) - 3.0*phi(3)/phi(20)) < 1e-12, "phi(3) reference failed");
rxDefault = Reaction(NominalRate=3.0);
assert(rxDefault.TemperatureResponse == "exponential", "default changed");
assert(abs(rxDefault.computeRate(20) - 3.0) < 1e-12, "exponential anchor failed");
fprintf("MATLAB CARDINAL TESTS PASS\n");
