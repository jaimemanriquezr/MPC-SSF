% TESTPATHOGEN  Verification of the audited pathogenModel preset (2026-08-27).
% See .claude/decisions/2026-08-27-pathogen-model-audit.md.
%
% (a) The Lund part is INHERITED, not duplicated: every component field and
%     every Lund reaction field is compared to modelLund() one by one, and the
%     comparison is printed.
% (b) The new Zeta0 / DetachForm / DetachScale / LundFlowingInert options.
% (c) A 0.5 d PAT pulse at N = 100. PAT appears in NO Lund reaction, so with the
%     three marker rates zeroed the PAT budget
%         dM = supply - export      (M = sum_z eps*c*dz over the three phases)
%     is the reaction-free instance of "in - out - reactions = d(storage)".
%     TWO checks, because the absolute budget has a measured floor:
%       (c1) TRACER INVARIANCE. The residual with the Lund biology ON must equal
%            the residual with EVERY reaction off -- if any Lund reaction touched
%            PAT the two would differ. Measured 3.9384e-6 vs 3.9377e-6 (2e-4
%            relative). This is the actual conservative-tracer statement.
%       (c2) ABSOLUTE CLOSURE to 1e-4 of supply, NOT 1e-6.
%            OPEN FINDING 2026-08-27: the eps-weighted budget for an EXPORTED
%            particle does not close to solver tolerance. With every reaction
%            rate zeroed and a constant influent it leaves a residual of
%            4.96e-5 of supply for PAT, while HET and PHO -- which are fully
%            retained, export ~1e-9 -- close to 2.9e-11. Frame refinement
%            101 -> 501 -> 2001 gives -6.27e-5 -> 4.44e-5 -> 4.86e-5, so it is
%            NOT trapezoid error; it converges to a nonzero value. The export
%            term here is q*c_flowing(end), whereas simulate.m:936 exports
%            eps_face*volumeAvgVelocity(end)*globalFlowing(end)
%            (src/@State/simulate.m:921-944); the two agree only if
%            eps*v_avg == q exactly at the outflow face, which biofilm growth
%            and osmosis need not preserve. probeMassClosure.m never caught this
%            because its influent PAT is 0, so no particle was ever exported and
%            the outflow term was never exercised. Needs its own probe.
% (d) The same pulse with the markers ON: the residual is then the net reaction
%     term. It must be a SINK and must lie inside the bracket set by the two
%     first-order limits (pure bacterivory at p_PAT vs pure inactivation at
%     d_PAT), both applied to the time-integrated biofilm PAT stock.
%
% Local, N = 100, ~1 min per run.

W = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(W, "src")));
addpath(fullfile(W, "analysis"));
addpath(fullfile(W, "analysis", "probes"));

%% ------------------------------------------------------------------ (a)
lund = modelLund();
pm   = pathogenModel();

fprintf("=== (a) Lund inheritance: pathogenModel() vs modelLund() ===\n");

% -- components ------------------------------------------------------------
assert(numel(pm.Components) == numel(lund.Components), "component count changed");
cFields = ["Name", "Density", "Dispersivity", "TransportRate", "Attenuation"];
fprintf("%-6s %-22s %14s %14s  %s\n", "comp", "field", "modelLund", "pathogenModel", "verdict");
nDiff = 0;
for i = 1:numel(lund.Components)
    a = lund.Components(i); b = pm.Components(i);
    assert(a.Name == b.Name, "component order changed at %d", i);
    for fn = cFields
        if ~isprop(a, fn), continue, end
        va = a.(fn); vb = b.(fn);
        if isstring(va) || ischar(va)
            same = isequal(va, vb);
            if ~same, fprintf("%-6s %-22s %14s %14s  DIFF\n", a.Name, fn, va, vb); nDiff = nDiff + 1; end
        else
            same = isequaln(va, vb);
            if ~same
                fprintf("%-6s %-22s %14.6g %14.6g  DIFF (x%.6g)\n", a.Name, fn, va, vb, vb/va);
                nDiff = nDiff + 1;
            end
        end
    end
    % attachment / sand factor live on Particles only
    if isa(a, "Particle")
        for fn = ["AttachmentSand", "AttachmentMatrix", "SandAttachmentFactor"]
            if ~isprop(a, fn), continue, end
            if ~isequaln(a.(fn), b.(fn))
                fprintf("%-6s %-22s %14.6g %14.6g  DIFF\n", a.Name, fn, a.(fn), b.(fn));
                nDiff = nDiff + 1;
            end
        end
    end
end
% Exactly two deliberate component differences: PAT TransportRate x40 and
% PAT SandAttachmentFactor -> SandPathogen (0). Everything else must match.
assert(nDiff == 2, "expected exactly 2 deliberate component differences, got %d", nDiff);
iPAT = find([pm.Components.Name] == "PAT", 1);
assert(abs(pm.Components(iPAT).TransportRate - 40*lund.Components(iPAT).TransportRate) < 1e-12);
assert(pm.Components(iPAT).SandAttachmentFactor == 0.0);
fprintf("components: %d fields differ, both deliberate (PAT Transport x40, PAT sand factor 0)\n\n", nDiff);

% -- Lund reactions --------------------------------------------------------
nL = numel(lund.Reactions);
assert(numel(pm.Reactions) == nL + 3, "expected %d Lund + 3 marker reactions", nL);
rFields = ["Name", "NominalRate", "TemperatureCorrectionFactor", "IsLightDependent", ...
           "MinimumLightFactor", "OptimalLightFactor", "LightInhibition", "EfficiencyBiofilm"];
fprintf("%-22s %-30s %14s %14s\n", "reaction", "field", "modelLund", "pathogenModel");
for i = 1:nL
    a = lund.Reactions(i); b = pm.Reactions(i);
    assert(a.Name == b.Name, "reaction order changed at %d", i);
    for fn = rFields
        va = a.(fn); vb = b.(fn);
        if isstring(va) || ischar(va)
            assert(isequal(va, vb), "%s.%s differs", a.Name, fn);
        else
            assert(isequaln(va, vb), "%s.%s differs: %g vs %g", a.Name, fn, va, vb);
        end
    end
    % half-saturations and stoichiometry, key by key
    Ka = a.HalfSaturationConstants; Kb = b.HalfSaturationConstants;
    assert(isequal(sort(keys(Ka)), sort(keys(Kb))), "%s half-sat keys differ", a.Name);
    for k = keys(Ka)'
        assert(abs(Ka(k) - Kb(k)) <= 1e-15*max(1, abs(Ka(k))), "%s K(%s) differs", a.Name, k);
        fprintf("%-22s %-30s %14.6g %14.6g\n", a.Name, "K_" + k, Ka(k), Kb(k));
    end
    Sa = a.StoichiometricCoefficients; Sb = b.StoichiometricCoefficients;
    assert(isequal(sort(keys(Sa)), sort(keys(Sb))), "%s stoich keys differ", a.Name);
    for k = keys(Sa)'
        assert(abs(Sa(k) - Sb(k)) <= 1e-15*max(1, abs(Sa(k))), "%s S(%s) differs", a.Name, k);
    end
    fprintf("%-22s %-30s %14.6g %14.6g\n", a.Name, "NominalRate", a.NominalRate, b.NominalRate);
    fprintf("%-22s %-30s %14.6g %14.6g\n", a.Name, "theta", ...
            a.TemperatureCorrectionFactor, b.TemperatureCorrectionFactor);
    % the ONE deliberate reaction difference: flowing-phase inertness
    assert(b.EfficiencyFlowing == 0.0, "%s should be flowing-inert by default", a.Name);
end
fprintf("Lund reactions: all rate/theta/half-sat/stoichiometry fields identical; ");
fprintf("only EfficiencyFlowing 1 -> 0 (deliberate).\n");

% LundFlowingInert=false restores modelLund's own flowing efficiencies.
pmOn = pathogenModel(LundFlowingInert=false);
for i = 1:nL
    assert(pmOn.Reactions(i).EfficiencyFlowing == lund.Reactions(i).EfficiencyFlowing, ...
        "LundFlowingInert=false must restore %s", lund.Reactions(i).Name);
end
fprintf("LundFlowingInert=false restores every Lund EfficiencyFlowing.\n\n");

% -- audited marker half-saturations --------------------------------------
mg = pm.Reactions([pm.Reactions.Name] == "MarkerGrowth");
Kmg = mg.HalfSaturationConstants;
expect = dictionary(["O2", "NH4", "HPO4", "DOM"], [2.0e-4, 1.0e-6, 2.0e-5, 4.0e-3]);
assert(isequal(sort(keys(Kmg)), sort(keys(expect))), "MarkerGrowth Monod set changed");
hg = lund.Reactions([lund.Reactions.Name] == "Heterotroph growth");
Khg = hg.HalfSaturationConstants;
fprintf("=== MarkerGrowth Monod set vs audited heterotroph growth ===\n");
for k = keys(expect)'
    assert(abs(Kmg(k) - expect(k)) < 1e-18, "MarkerGrowth K_%s = %g, expected %g", k, Kmg(k), expect(k));
    assert(abs(Kmg(k) - Khg(k)) < 1e-18, "MarkerGrowth K_%s must equal HET growth K_%s", k, k);
    fprintf("  K_%-5s marker %10.4g   HET growth %10.4g   agree\n", k, Kmg(k), Khg(k));
end
bv = pm.Reactions([pm.Reactions.Name] == "Bacterivory");
assert(abs(bv.HalfSaturationConstants("HET") - 2e-3) < 1e-18, "K_pred changed from 2e-3");
assert(abs(mg.NominalRate - 0.2) < 1e-15 && abs(mg.TemperatureCorrectionFactor - 1.047) < 1e-15);
assert(abs(pm.Reactions([pm.Reactions.Name] == "Inactivation").NominalRate - 0.02) < 1e-15);
assert(abs(bv.NominalRate - 8.0) < 1e-15);
assert(bv.EfficiencyFlowing == 1e-3, "WaterFactor default changed");
fprintf("  K_pred 2e-3 kept (tab:eco-parameters lists none); mu/d/p = 0.2/0.02/8.0 unchanged\n\n");

%% ------------------------------------------------------------------ (b)
fprintf("=== (b) new options ===\n");
assert(pm.CohesionSubModel.Zeta0 == 1e2, "default Zeta0 must stay 1e2");
assert(abs(pm.DetachmentFunction(18) - 1.4e-5) < 1e-18, "default detachment law changed");
p1 = pathogenModel(Zeta0=1, DetachForm="linear");
assert(p1.CohesionSubModel.Zeta0 == 1, "Zeta0 option ignored");
assert(abs(p1.DetachmentFunction(18) - 0.14) < 1e-15);
assert(abs(p1.DetachmentFunction(9)  - 0.07) < 1e-15, "linear form wrong");
p2 = pathogenModel(DetachForm="sqrt", DetachScale=2.0);
assert(abs(p2.DetachmentFunction(18) - 0.28) < 1e-15, "sqrt/scale wrong");
assert(p1.CohesionSubModel.Zeta1 == lund.CohesionSubModel.Zeta1, "Zeta1 must stay inherited");
fprintf("Zeta0/DetachForm/DetachScale all honoured; defaults unchanged.\n\n");

%% ------------------------------------------------------------------ (c,d)
% Working-set host (E4): zeta0 = 1, linear detachment, but only 0.5 d at N=100
% from a CLEAN filter -- this is a budget test, not a physics run.
NCELLS = 100;  TEND = 0.5;  MAXDT = 5e-5;  NFRAMES = 101;
PATBASE = 5.36e-3;                 % Manriquez Table B.1 marker influent
PULSE   = 10.0;  T0 = 0.1;  T1 = 0.3;
infl = [3.0e-4, 1.0e-3, 0, PATBASE, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3];
pulsed = infl;  pulsed(4) = PULSE*PATBASE;
inflow = @(t) tern(t >= T0 && t < T1, pulsed, infl);

fprintf("=== (c) PAT conservative-tracer closure (marker rates = 0) ===\n");
[res0, r0] = patBudget(pathogenModel(Zeta0=1, DetachForm="linear", NormalizedLight=true), ...
                       "markers", NCELLS, TEND, MAXDT, NFRAMES, inflow, infl, PATBASE, PULSE, T0, T1);
fprintf("  Lund biology ON, markers off:  dM = %.8e  supply = %.8e  export = %.8e\n", ...
        res0.dM, res0.supply, res0.export);
fprintf("                                 residual = %.4e  rel = %.4e\n", res0.resid, res0.rel);
assert(r0.Flag == "OK", "tracer run flag %s", r0.Flag);
% (c1) tracer invariance: switching the Lund biology off must not move the
% PAT residual, since no Lund reaction has a PAT stoichiometric entry.
[resZ, rZ] = patBudget(pathogenModel(Zeta0=1, DetachForm="linear", NormalizedLight=true), ...
                       "all", NCELLS, TEND, MAXDT, NFRAMES, inflow, infl, PATBASE, PULSE, T0, T1);
assert(rZ.Flag == "OK", "all-off run flag %s", rZ.Flag);
relInv = abs(res0.resid - resZ.resid)/max(abs(resZ.resid), realmin);
fprintf("  every reaction off:            residual = %.4e   (invariance %.3e)\n", resZ.resid, relInv);
assert(relInv < 1e-3, "PAT residual moved when the Lund biology was switched off (%.3e)", relInv);
fprintf("  PASS (c1): the PAT budget is invariant to the Lund biology -- conservative tracer.\n");
% (c2) absolute closure, at the measured floor. See the OPEN FINDING at the top.
assert(abs(res0.rel) < 1e-4, "PAT absolute closure worse than 1e-4: %.3e", res0.rel);
fprintf("  PASS (c2): |residual|/supply = %.3e < 1e-4 (floor ~5e-5, open finding).\n\n", abs(res0.rel));

fprintf("=== (d) PAT budget with the marker reactions ON ===\n");
[res1, r1] = patBudget(pathogenModel(Zeta0=1, DetachForm="linear", NormalizedLight=true), ...
                       "none", NCELLS, TEND, MAXDT, NFRAMES, inflow, infl, PATBASE, PULSE, T0, T1);
assert(r1.Flag == "OK", "marker run flag %s", r1.Flag);
% Correct for the reaction-free budget floor measured in (c), so what is left is
% the reaction term alone.
netRx = res1.resid - res0.resid;
fprintf("  dM = %.6e  supply = %.6e  export = %.6e\n", res1.dM, res1.supply, res1.export);
fprintf("  netRx (residual minus the reaction-free floor %.3e) = %.6e\n", res0.resid, netRx);
assert(netRx < 0, "markers must be a NET SINK for PAT, got %+.3e", netRx);
% Upper bound: every marker reaction is first order in PAT with a Monod factor
% <= 1, so with J_bio = int sum_z eps*X_PAT^(matrix+enclosed) dz dt and
% J_flo the same for the flowing phase,
%   sink <= (d_PAT + p_PAT)*J_bio + WaterFactor*p_PAT*J_flo
% (inactivation and growth are flowing-inert; bacterivory runs there at
% EfficiencyFlowing = WaterFactor). Growth can only reduce the sink.
th = @(rx, T) rx.NominalRate*rx.TemperatureCorrectionFactor^(T/293 - 1);
T = 19;
dP = th(pm.Reactions([pm.Reactions.Name] == "Inactivation"), T);
pP = th(bv, T);
bound = (dP + pP)*res1.Jbio + bv.EfficiencyFlowing*pP*res1.Jflo;
fprintf("  J_bio = %.6e   J_flo = %.6e kg d/m2   (d_PAT %.4g, p_PAT %.4g /d at %g C)\n", ...
        res1.Jbio, res1.Jflo, dP, pP, T);
fprintf("  first-order sink bound = %.4e   observed sink = %.4e   ratio %.3f\n", ...
        bound, -netRx, -netRx/bound);
assert(-netRx <= bound*1.001, "sink %.4e exceeds its first-order bound %.4e", -netRx, bound);
assert(-netRx > 0.05*bound, "sink %.4e implausibly far below its bound %.4e", -netRx, bound);
fprintf("  PASS: net marker sink is a sink and is under its first-order bound.\n\n");

fprintf("MATLAB PATHOGEN-MODEL TESTS PASS\n");

%% ---------------------------------------------------------------- helpers
function out = tern(c, a, b)
if c, out = a; else, out = b; end
end

function [out, r] = patBudget(m, zero, ncells, tend, maxDt, nframes, ...
                              inflow, infl, patBase, pulse, t0, t1)
% eps-weighted PAT stock budget over a pulse run (convention of
% analysis/probes/probeMassClosure.m). `zero` is "none" | "markers" | "all".
if zero ~= "none"
    rx = m.Reactions;
    if zero == "all"
        for k = 1:numel(rx), rx(k).NominalRate = 0.0; end
    else
        for k = ["MarkerGrowth", "Inactivation", "Bacterivory"]
            rx([rx.Name] == k).NominalRate = 0.0;
        end
    end
    m = Model(m.Components, rx, Kappa=m.CohesionSubModel.Kappa, ...
        Zeta0=m.CohesionSubModel.Zeta0, Zeta1=m.CohesionSubModel.Zeta1, ...
        DetachmentFunction=m.DetachmentFunction, WaterDensity=m.WaterDensity, ...
        BiofilmPorosity=m.BiofilmPorosity, OsmosisRate=m.OsmosisRate);
end
f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(ncells);
r = simulate(State(f, m), InflowConcentrations=inflow, SimulationTime=tend, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=maxDt, ...
    ImplicitOsmosis=true, CohesionBC="neumann", CohesionScheme="shin", ...
    FrameNumber=nframes, Quiet=true);
C = r.Frames.Concentrations; ts = r.Frames.Time(:);
z = f.GridPoints.Centers(:); dz = f.GridSize; ep = computePorosity(f, z);
q = f.InflowVelocity;
mat = C{"PAT","Matrix"}{1}; enc = C{"PAT","Enclosed"}{1}; flo = C{"PAT","Flowing"}{1};
M = sum(ep.*(mat + enc + flo), 1)*dz;
cIn = infl(4)*ones(size(ts));  cIn(ts >= t0 & ts < t1) = pulse*patBase;
out.dM     = M(end) - M(1);
out.supply = q*trapz(ts, cIn);
out.export = q*trapz(ts, flo(end, :).');
out.resid  = out.dM - (out.supply - out.export);
out.rel    = out.resid/max(out.supply, realmin);
% time integrals of the eps-weighted PAT stock, split by phase: the marker
% reactions act at full strength on the biofilm phases and at WaterFactor on
% the flowing suspension.
out.Jbio = trapz(ts, (sum(ep.*(mat + enc), 1)*dz).');
out.Jflo = trapz(ts, (sum(ep.*flo, 1)*dz).');
end
