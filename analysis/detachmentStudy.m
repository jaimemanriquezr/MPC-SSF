function outdir = detachmentStudy(options)
% DETACHMENTSTUDY  Detachment prefactor as an uncertain parameter (Jaime,
% 2026-08-20: "detachment is a simulation parameter worth exploring").
%
% Context: the legacy repo used b_det = 0.14*sqrt(v/18) for biofilm runs but
% 1.4e-5*sqrt(v/18) for pathogen runs; under the corrected model the latter
% clogs within the 10-day marker window. The 1.5-day OAT campaign found
% `detach_scale` inert -- this study shows its long-horizon role: sweep the
% prefactor over four orders of magnitude, 10-day constant marker feed from
% the mature 30-day summer state (manuscript cache), and record log removal,
% peak biofilm fraction, and clogging.
arguments
    options.NCells (1,1) {mustBeNumeric} = 100;
    options.MaxDt (1,1) {mustBeNumeric} = 3e-6;
    options.PatLen (1,1) {mustBeNumeric} = 10;
end
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, "results", "manuscript", "X2_detachment");
if ~isfolder(outdir), mkdir(outdir); end
cache = fullfile(here, "results", "manuscript", "cache", ...
    sprintf("mature30_summer_N%d.mat", options.NCells));
assert(isfile(cache), "mature cache missing; run manuscriptExperiments(""E3"") first");
st = load(cache).st;

prefactors = [1.4e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.14, 0.5];
patNom = 5.3616e-3;
infl = [2.68e-3, 1.00e-2, 0.0, patNom, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
header = ["prefactor", "flag", "t_final", "logremoval_final", "logremoval_mean", "phib_peak"];
rows = cell(length(prefactors), length(header));
for i = 1:length(prefactors)
    b0 = prefactors(i);
    m = pathogenModel(PhototrophRespiration=0.55, PGExcess=true);
    m = Model(m.Components, m.Reactions, Kappa=1e-7, ...
        Zeta0=m.CohesionSubModel.Zeta0, Zeta1=m.CohesionSubModel.Zeta1, ...
        DetachmentFunction=@(v) b0*sqrt(abs(v)/18), WaterDensity=m.WaterDensity, ...
        BiofilmPorosity=m.BiofilmPorosity, OsmosisRate=m.OsmosisRate);
    f = SandFilter(Temperature=19);
    f = f.addGridPoints(options.NCells);
    f.LightIrradiation = @(t) max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);
    s = State(f, m);
    s.GlobalConcentration.Matrix = st.Matrix;
    s.GlobalConcentration.EnclosedParticles = st.EnclosedParticles;
    s.GlobalConcentration.FlowingParticles = st.FlowingParticles;
    s.GlobalConcentration.EnclosedLiquids = st.EnclosedLiquids;
    s.GlobalConcentration.FlowingLiquids = st.FlowingLiquids;
    s.EnclosedWaterVolume = st.EnclosedWaterVolume;
    s.Velocity.Biofilm = st.VelocityBiofilm;
    s.Time = 30;
    r = simulate(s, InflowConcentrations=infl, SimulationTime=options.PatLen, ...
        TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=options.MaxDt, ...
        FrameNumber=48, ImplicitOsmosis=true, Quiet=true);
    C = r.Frames.Concentrations;
    ts = r.Frames.Time(:);
    patOut = C{"PAT", "Flowing"}{1}(end, :)';
    F = log10(patNom./max(patOut, 1e-30));
    densityL = mean([m.Liquids.Density]); densityP = mean([m.Particles.Density]);
    phiB = C{"Water","Enclosed"}{1}/densityL;
    for nm = [m.Particles.Name], phiB = phiB + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/densityP; end
    for nm = [m.Liquids.Name], phiB = phiB + C{nm,"Enclosed"}{1}/densityL; end
    rows(i, :) = {b0, string(r.Flag), r.TimeFinal, F(end), mean(F(ts > 30.5)), max(phiB, [], "all")};
    writematrix([ts, patOut, F], fullfile(outdir, sprintf("logremoval_b%g.csv", b0)));
    fprintf("  b0=%-8g flag=%-8s F_end=%.3f F_mean=%.3f phib=%.3f\n", ...
        b0, r.Flag, F(end), mean(F(ts > 30.5)), max(phiB, [], "all"));
end
writetable(cell2table(rows, 'VariableNames', cellstr(header)), fullfile(outdir, "summary.csv"));
fprintf("wrote %s\n", fullfile(outdir, "summary.csv"));
end
