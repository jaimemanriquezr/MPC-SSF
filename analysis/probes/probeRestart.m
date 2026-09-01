function out = probeRestart(opts)
% PROBERESTART  Does chaining a run reproduce a single continuous run?
%
% GATE for every long job. Runs >10 d are to be executed 10 d at a time, feeding
% the final state of one leg into the next. That is only legitimate if the
% snapshot carries the COMPLETE state. State holds:
%   GlobalConcentration.{Matrix, EnclosedParticles, FlowingParticles,
%                        EnclosedLiquids, FlowingLiquids}
%   EnclosedWaterVolume, Velocity.{Biofilm, Flowing}, Time
% stateFromFrame/rehomeState restore all of these except Velocity.Flowing, which
% simulate recomputes from velBiofilm and phiBiofilm before the loop.
%
% The risk is silent: stateFromFrame's `grab` substitutes ZEROS for any component
% stored as a scalar, and State's constructor zero-fills everything it is not
% given. Either would produce a restart that runs happily from partly-empty
% initial data. So this compares a chained run against a continuous one and
% reports the difference field by field -- an assertion, not an inspection.
arguments
    opts.NCells (1,1) double = 100
    opts.Days (1,1) double = 6
    opts.Legs (1,1) double = 2
    opts.MaxDt (1,1) double = 3e-6
    opts.Kappa (1,1) double = 1e-6
    opts.Zeta1 (1,1) double = 0.27
    opts.Tol (1,1) double = 1e-3
    % Continuation legs restart the adaptive stepper. At 1e-8 (the cold-start
    % value) each leg re-ramps from scratch, which perturbs the dt sequence at the
    % leg boundary; seeding near the working dt removes most of that transient.
    opts.ContinuationDt (1,1) double = 0
end
here = fileparts(mfilename("fullpath")); W = fileparts(fileparts(here));
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis")); addpath(here);
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=opts.Zeta1, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
common = {"InflowConcentrations", infl, "TimeStep", "adaptive", ...
          "AdaptiveInitialDt", 1e-8, "AdaptiveMaxDt", opts.MaxDt, ...
          "ImplicitOsmosis", true, "ImplicitDispersion", true, "Quiet", true};

fprintf("restart gate: N=%d, %g d as 1 leg vs %d legs, zeta_1=%g, MaxDt=%g\n", ...
    opts.NCells, opts.Days, opts.Legs, opts.Zeta1, opts.MaxDt);

% --- reference: one continuous run -----------------------------------------
t0 = tic;
rRef = simulate(State(f, m), common{:}, SimulationTime=opts.Days, FrameNumber=4);
checkRunFlag(rRef, "continuous reference");
fprintf("  continuous: %s, %.0f s\n", rRef.Flag, toc(t0));

% --- chained: Legs legs of Days/Legs each ----------------------------------
leg = opts.Days/opts.Legs;  s = State(f, m);  t = 0;  zeroed = string.empty;
t1 = tic;
for L = 1:opts.Legs
    dt0 = 1e-8;
    if L > 1 && opts.ContinuationDt > 0, dt0 = opts.ContinuationDt; end
    rL = simulate(s, common{:}, SimulationTime=leg, FrameNumber=4, AdaptiveInitialDt=dt0);
    checkRunFlag(rL, sprintf("chained leg %d/%d", L, opts.Legs));
    st = stateFromFrame(rL, numel(rL.Frames.Time));
    if ~isempty(st.allZero), zeroed = [zeroed, st.allZero]; end %#ok<AGROW>
    t = t + leg;  s = rehomeState(f, m, st, t);
    fprintf("  leg %d/%d done (t = %g d), all-zero fields: %s\n", L, opts.Legs, t, ...
        ternaryLocal(isempty(st.allZero), "none", strjoin(st.allZero, ", ")));
end
rChain = simulate(s, common{:}, SimulationTime=0, FrameNumber=2); %#ok<NASGU>
fprintf("  chained: %.0f s\n", toc(t1));

% --- compare the two final states field by field ---------------------------
A = stateFromFrame(rRef, numel(rRef.Frames.Time));
B = st;
out.fields = ["Matrix","EnclosedParticles","FlowingParticles","EnclosedLiquids", ...
              "FlowingLiquids","EnclosedWaterVolume","VelocityBiofilm"];
out.relDiff = zeros(1, numel(out.fields));
fprintf("\n  %-20s %12s %12s %12s\n", "field", "max|ref|", "max|chain|", "max rel diff");
for q = 1:numel(out.fields)
    a = A.(out.fields(q)); b = B.(out.fields(q));
    d = max(abs(a(:) - b(:)))/max(max(abs(a(:))), realmin);
    out.relDiff(q) = d;
    fprintf("  %-20s %12.4e %12.4e %12.3e%s\n", out.fields(q), max(abs(a(:))), max(abs(b(:))), d, ...
        ternaryLocal(d > opts.Tol, "   <-- FAIL", ""));
end
out.zeroedFields = unique(zeroed);
out.worst = max(out.relDiff);
out.pass = out.worst <= opts.Tol && isempty(out.zeroedFields);
fprintf("\n  worst field difference: %.3e (tol %.1e)\n", out.worst, opts.Tol);
fprintf("  fields zeroed by a snapshot: %s\n", ternaryLocal(isempty(out.zeroedFields), "none", strjoin(out.zeroedFields, ", ")));
fprintf("  GATE: %s\n", ternaryLocal(out.pass, "PASS -- chaining is safe", "FAIL -- do NOT chain"));
save(fullfile(S, sprintf("restart_gate_n%d_%gd_%dlegs.mat", opts.NCells, opts.Days, opts.Legs)), "out", "-v7");
end

function v = ternaryLocal(c, a, b)
if c, v = a; else, v = b; end
end
