function s = seedState(filter, model, phiBed)
% SEEDSTATE  A State carrying a DEVELOPED biofilm, built analytically.
%
%   s = seedState(filter, model)           % phi_b = 0.25 at the bed surface
%   s = seedState(filter, model, 0.4)
%
% WHY THIS EXISTS. The regression configs run for 0.05 d. Started from the zero
% State, that leaves the biofilm matrix at ~1e-9 and the enclosed phases at
% ~1e-18, with half the domain still exactly zero -- so the cohesion solve, the
% detachment law, the attachment split and every Monod term are exercised on
% numerical nothing. Measured, not assumed: the first capture of
% baseline_theta_mixed had max|HET/Enclosed| = 5.2e-19. A gate like that cannot
% see a bug in Solver A, which is precisely the code the cleanup moves.
%
% Growing a real biofilm instead would cost days of simulated time per config.
% Seeding one costs nothing and is exactly reproducible, which is all the gate
% needs -- this profile is a plausible state, not a physical prediction, and
% nothing but the gate should ever read it.
%
% CONSTRUCTION. A target phi_b profile is split the way the osmosis relaxation
% would leave it at equilibrium (rhsEnclosedWater = (beta*phi_b - phi_e)/tau, so
% phi_e -> beta*phi_b):
%
%   phi_matrix = (1 - beta) * phi_b          the skeleton
%   phi_enc    = beta * phi_b                enclosed water + enclosed solutes
%   phi_f      = 1 - phi_b                   flowing
%
% so the state starts near the osmotic fixed point rather than relaxing violently
% through the first few steps. Solutes AND suspended particles carry the influent
% concentrations in every phase, which likewise avoids a startup transient -- see
% the note above `pconc` for what a load above the influent did to the gate.

arguments
    filter  mpcssf.SandFilter
    model   mpcssf.Model
    phiBed (1,1) double = 0.25
end

z  = filter.GridPoints.Centers(:);
N  = numel(z);
n0 = filter.GridZero;

kP = numel(model.Particles);
kL = numel(model.Liquids);
rhoP = mean([model.Particles.Density]);
rhoL = mean([model.Liquids.Density]);
beta = model.BiofilmPorosity;

% Target phi_b: zero in the supernatant (cells 1..n0, above z = 0), then an
% exponential decay with depth into the bed. The Schmutzdecke shape -- most of
% the biomass in the top centimetre or two -- rather than a uniform slab, so the
% cohesive flux sees a real gradient to act on.
phib = zeros(N, 1);
bed = (n0+1):N;
depth = z(bed) - z(n0+1);
phib(bed) = phiBed * exp(-abs(depth)/0.02);

% Species split of the matrix. Fixed, arbitrary, and documented as such: the
% gate needs a repeatable state, not a calibrated one. PAT is kept small because
% it is the manuscript's headline output and a large seed would mask a change in
% it; PG (present only in the "wolf" respiration branch) starts empty, since its
% whole point is to be filled by photosynthesis.
frac = zeros(1, kP);
pnames = [model.Particles.Name];
for k = 1:kP
    switch pnames(k)
        case "HET",  frac(k) = 0.50;
        case "PHO",  frac(k) = 0.30;
        case "POM",  frac(k) = 0.19;
        case "PAT",  frac(k) = 0.01;
        otherwise,   frac(k) = 0.0;    % PG and anything a future preset adds
    end
end
frac = frac/sum(frac);

% SUSPENDED particle concentrations, per unit volume of the phase holding them.
% These are the influent values the configs feed, listed here for the same reason
% lconc is: the seed is a plausible state, and a suspended load far above the
% influent is not one.
%
% Until 2026-09-07 this was `1e-3 * rhoP * frac` (flowing) and `1e-2 * rhoP *
% frac` (enclosed), which gave 0.559 kg/m^3 of flowing HET against an influent of
% 3.0e-4 -- about 1900x -- and 0.212 kg/m^3 of POM out of nothing. That was inert
% at AttachmentMatrix = 5.47e2, but at 2.735e4 (since 2026-09-05) the load
% plates onto the bed-surface cell within hours and trips the clogging guard: 15
% of the 18 configs aborted CLOGGED at t ~ 0.037-0.039 d, and refining 30 -> 100
% cells moved that to 0.0148 d, because a thinner surface cell concentrates the
% same captured mass further. Measured 2026-09-07; with the load at the influent
% values the same seeded film runs the full 0.05 d at Flag = OK. The clogging was
% a defect in THIS FILE, not evidence about the attachment rate.
%
% POM has no influent entry, so it takes an influent-SCALE value instead of zero,
% which keeps the flowing- and enclosed-POM terms doing arithmetic. PG stays at 0
% for the reason its matrix fraction does: photosynthesis is what fills it.
pconc = zeros(1, kP);
for k = 1:kP
    switch pnames(k)
        case "HET",  pconc(k) = 3.00e-4;
        case "PHO",  pconc(k) = 1.00e-3;
        case "POM",  pconc(k) = 1.00e-3;   % no influent entry; influent-scale
        case "PAT",  pconc(k) = 5.36e-3;
        otherwise,   pconc(k) = 0.0;       % PG and anything a future preset adds
    end
end

% Solute concentrations, per unit volume of the phase holding them.
lnames = [model.Liquids.Name];
lconc = zeros(1, kL);
for k = 1:kL
    switch lnames(k)
        case "O2",   lconc(k) = 9.10e-3;
        case "IC",   lconc(k) = 1.70e-3;
        case "NH4",  lconc(k) = 2.00e-5;
        case "HPO4", lconc(k) = 9.30e-6;
        case "DOM",  lconc(k) = 1.00e-3;
        otherwise,   lconc(k) = 0.0;
    end
end

phiMatrix = (1 - beta)*phib;
phiEnc    = beta*phib;
phiFlow   = 1 - phib;

s = mpcssf.State(filter, model);

s.GlobalConcentration.Matrix = rhoP * phiMatrix .* frac;

% Enclosed particles: the suspended load inside the biofilm. Nonzero is the point
% -- it is what makes the enclosed attachment term attE and the particle transfer
% term transP do arithmetic.
s.GlobalConcentration.EnclosedParticles = phiEnc .* pconc;
s.GlobalConcentration.EnclosedLiquids   = phiEnc .* lconc;

% Enclosed water takes whatever the solutes and suspended particles leave, so
% that phi_e comes out at beta*phi_b exactly.
s.EnclosedWaterVolume = phiEnc ...
    - sum(s.GlobalConcentration.EnclosedParticles, 2)/rhoP ...
    - sum(s.GlobalConcentration.EnclosedLiquids, 2)/rhoL;

s.GlobalConcentration.FlowingParticles = phiFlow .* pconc;
s.GlobalConcentration.FlowingLiquids   = phiFlow .* lconc;

% Fail loudly rather than hand the solver a state it will reject on step one.
phibCheck = s.EnclosedWaterVolume ...
    + sum(s.GlobalConcentration.Matrix, 2)/rhoP ...
    + sum(s.GlobalConcentration.EnclosedParticles, 2)/rhoP ...
    + sum(s.GlobalConcentration.EnclosedLiquids, 2)/rhoL;
if max(abs(phibCheck - phib)) > 1e-12
    error("seedState:inconsistent", ...
        "seeded phi_b differs from target by %.3e", max(abs(phibCheck - phib)));
end
if any(s.EnclosedWaterVolume < 0)
    error("seedState:negativeWater", "seeded enclosed water is negative");
end
if max(phib) >= 0.99
    error("seedState:clogged", "seeded phi_b = %.3f would trip the clogging guard", max(phib));
end
end
