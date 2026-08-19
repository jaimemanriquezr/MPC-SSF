function outdir = phoBloomStudy(options)
% PHOBLOOMSTUDY  Effect of phototroph inflow concentration and growth/death
% rates on the oxygen field (algae-bloom hypothesis, 2026-08-18).
%
% Hypothesis (Jaime): the simulated algae bloom is too strong and produces the
% unrealistic oxygen increase flagged by Reviewer 1. This study grows a filter
% from a clean state for TSim days (the regime where the bloom and the O2 rise
% develop) and sweeps, one at a time:
%   * influent_PHO : influent phototroph concentration, nominal 1.00e-2 kg/m3
%   * mu_PHO       : phototroph growth rate, nominal 5.5 /d
%   * d_PHO        : phototroph death rate, nominal 0.4 /d
% plus two combined "suppressed bloom" runs. QoIs per run (final frame unless
% noted): max O2 anywhere over the whole run / influent O2 (supersaturation
% ratio), outflow O2, total phototroph biofilm mass, peak biofilm fraction, and
% the full final O2/PHO profiles (per-run CSV).
arguments
    options.TSim (1,1) {mustBeNumeric} = 3.0;
    options.NCells (1,1) {mustBeNumeric} = 30;
    options.NFrames (1,1) {mustBeNumeric} = 60;
    options.MaxDt (1,1) {mustBeNumeric} = 3e-6;
    % Override the phototroph dark-respiration light floor (MinimumLightFactor
    % on "Phototroph growth"; modelLund default 0.01). 0 = photosynthesis fully
    % off in darkness. NaN = leave the preset untouched.
    options.DarkRespiration (1,1) {mustBeNumeric} = NaN;
    % > 0 enables the phototroph metabolic split (photosynthesis + maintenance
    % respiration; see modelLund.m). Campos2006 kra avg = 0.276 /d.
    options.Respiration (1,1) {mustBeNumeric} = 0.0;
    options.OutDirName (1,1) string = "pho_bloom";
end

influentBase = [2.68e-3, 1.00e-2, 0.0, 5.36e-3, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
o2In = influentBase(5);
muNom = 5.5;  dNom = 0.4;  phoInNom = influentBase(2);

% One-at-a-time sweeps + combined suppressed-bloom runs. factor applies to the
% named parameter; nominal (factor 1) is shared as the baseline run.
runs = struct("label", {}, "phoIn", {}, "mu", {}, "d", {});
add = @(label, phoIn, mu, d) struct("label", string(label), "phoIn", phoIn, "mu", mu, "d", d);
runs(end+1) = add("baseline", phoInNom, muNom, dNom);
for fac = [0.1, 0.3, 3.0]
    runs(end+1) = add(sprintf("influent_PHO_x%g", fac), fac*phoInNom, muNom, dNom);
end
for fac = [0.25, 0.5, 2.0]
    runs(end+1) = add(sprintf("mu_PHO_x%g", fac), phoInNom, fac*muNom, dNom);
end
for fac = [0.5, 2.0, 4.0]
    runs(end+1) = add(sprintf("d_PHO_x%g", fac), phoInNom, muNom, fac*dNom);
end
runs(end+1) = add("suppressed_mild",   0.3*phoInNom, 0.5*muNom, 2.0*dNom);
runs(end+1) = add("suppressed_strong", 0.1*phoInNom, 0.25*muNom, 4.0*dNom);

here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, "results", options.OutDirName);
if ~isfolder(outdir), mkdir(outdir); end
fprintf("PHO bloom study: %d runs, tsim=%g ncells=%d\n", length(runs), options.TSim, options.NCells);

header = ["run", "influent_PHO", "mu_PHO", "d_PHO", "flag", "o2_max_ratio", ...
          "o2_out_final", "o2_min_final", "pho_biofilm_mass", "phib_peak"];
rows = cell(length(runs), length(header));
for i = 1:length(runs)
    q = runs(i);
    f = SandFilter();
    f = f.addGridPoints(options.NCells);
    m = pathogenModel(PhototrophRespiration=options.Respiration);
    rx = m.Reactions;
    for j = 1:length(rx)
        if rx(j).Name == "Phototroph growth", rx(j).NominalRate = q.mu; end
        if rx(j).Name == "Phototroph death",  rx(j).NominalRate = q.d; end
        if ~isnan(options.DarkRespiration) && rx(j).Name == "Phototroph growth"
            rx(j).MinimumLightFactor = options.DarkRespiration;
        end
    end
    m.Reactions = rx;
    infl = influentBase;  infl(2) = q.phoIn;
    if options.Respiration > 0
        infl = [infl(1:4), 0.0, infl(5:end)];   % PG = 0 in the influent
    end

    r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=options.TSim, ...
        TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=options.MaxDt, ...
        FrameNumber=options.NFrames, ImplicitOsmosis=true, Quiet=true);

    C = r.Frames.Concentrations;
    o2Flow = C{"O2", "Flowing"}{1};                      % N x frames
    o2EnclConc = C{"O2", "Enclosed"}{1};                 % biofilm-phase O2 (global)
    phoBiofilm = C{"PHO", "Matrix"}{1} + C{"PHO", "Enclosed"}{1};
    densityL = mean([m.Liquids.Density]);
    densityP = mean([m.Particles.Density]);
    phiW = C{"Water", "Enclosed"}{1}/densityL;
    pNames = [m.Particles.Name];  lNames = [m.Liquids.Name];
    phiB = phiW;
    for nm = pNames
        phiB = phiB + (C{nm, "Matrix"}{1} + C{nm, "Enclosed"}{1})/densityP;
    end
    for nm = lNames
        phiB = phiB + C{nm, "Enclosed"}{1}/densityL;
    end

    dz = f.GridSize;
    o2MaxRatio = max(o2Flow, [], "all")/o2In;
    o2OutFinal = o2Flow(end, end);
    o2MinFinal = min(o2Flow(:, end));
    phoMass = sum(phoBiofilm(:, end))*dz;
    phibPeak = max(phiB(:, end));

    rows(i, :) = {q.label, q.phoIn, q.mu, q.d, string(r.Flag), o2MaxRatio, ...
                  o2OutFinal, o2MinFinal, phoMass, phibPeak};
    writematrix([f.GridPoints.Centers(:), o2Flow(:, end), o2EnclConc(:, end), ...
                 phoBiofilm(:, end), phiB(:, end)], ...
                fullfile(outdir, "profiles_" + q.label + ".csv"));
    fprintf("  %-22s flag=%-8s O2max/O2in=%.3f  O2out=%.4g  PHOmass=%.4g  phib=%.4f\n", ...
            q.label, r.Flag, o2MaxRatio, o2OutFinal, phoMass, phibPeak);
end

writetable(cell2table(rows, 'VariableNames', cellstr(header)), ...
           fullfile(outdir, "summary.csv"));
fprintf("wrote %s\n", fullfile(outdir, "summary.csv"));
end
