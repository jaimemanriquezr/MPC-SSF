# MATLAB Log-OAT campaign + PHO bloom study

## Goal

Run the decided-methodology sensitivity campaign **in MATLAB** (MATLAB anchors as golden
reference; Julia already ran the mirrored design), with two design changes Jaime ordered
on 2026-08-18:
1. `velocity` excluded (fixed operating condition, see
   `.claude/decisions/2026-08-18-velocity-fixed-operating-condition.md`) → 26 parameters.
2. `beta_porosity` no longer gap ×2/×½ (old scheme left the [0.95, 0.99] interval at
   β=0.995): perturbed runs at **β = 0.98 and β = 0.95** (gap 0.02, 0.05), log-secant
   slope over ln(gap), denominator ln(2.5). Baseline stays β = 0.99.

Plus a **new PHO bloom study**: influent PHO concentration, `mu_PHO`, `d_PHO` swept with
O₂ QoIs, to test Jaime's suspicion that the simulated algae bloom is too strong and
causes the unrealistic O₂ increase (Reviewer 1 point 3).

Done = three `measures.csv` (startup/pulse/flowstep) + PHO study CSV, produced by MATLAB,
parameters documented, results comparable to `julia/analysis/results/log_oat_*`.

## Where

Worktree `.claude/worktrees/agent-a5f45f8cf63f54f77` (branch `matlab-claude` — the only
MATLAB tree with the pathogen port, light floor, Celsius temperature, adaptive CFL).
`julia-port`'s MATLAB tree lacks all four (0 src/ commits ahead).

## Steps

1. **Port `implicit_osmosis` to `src/@State/simulate.m`** (from `implicit-osmosis.md` /
   `julia/src/simulate.jl:403-518`): new `ImplicitOsmosis` option; damp
   `rhsEnclosedWater` by `1/(1+dt·kosm)`, `kosm=(1−β)/τ`, at both use sites (Solver A
   volume RHS, φW update); drop `1/tau` from CFL `w_b(3)` when implicit.
   → verify: smoke run reaches dt = AdaptiveMaxDt = 3e-6 instead of ~1e-7.
2. **`src/presets/pathogenModel.m`** — MATLAB port of `pathogen_repro.jl:pathogen_model`:
   modelLund + PAT Transport ×40 + SandAttachmentFactor + 3 marker reactions
   (MarkerGrowth 0.2/θ1.047, Inactivation 0.02/θ1.08, Bacterivory 8.0/θ1.08 K_HET=2e-3,
   EfficiencyFlowing=water_factor) + Lund reactions inert in flowing + Zeta0=1e2 +
   detachment 1.4e-5·√(|v|/18). → verify: component/reaction counts, stoichiometry spot
   check vs Julia.
3. **`analysis/logOatCampaign.m`** — mirror of `log_oat_sensitivity.jl` `main()`:
   26 params, scenarios startup/pulse/flowstep (pulse 10× over [0.1,0.3), surge ×1.4,
   tmature=3.0, tpost=1.5, ncells=30, nframes=60), measures I_rms/I_max/D_min/I_min/asym,
   per-param curves CSV. β special-cased. → verify: smoke run (tiny times) end-to-end;
   orientation assert `centers(end)>0`.
4. **`analysis/phoBloomStudy.m`** — clean start, tsim=3.0 (the regime where the bloom
   and O₂ rise develop), sweep influent_PHO ∈ nominal×{0.1,0.3,1,3}, mu_PHO ∈
   5.5×{0.25,0.5,1,2}, d_PHO ∈ 0.4×{0.5,1,2,4}, plus 2 combined "suppressed bloom" runs.
   QoIs: max O₂ anywhere / influent O₂ (supersaturation ratio), outflow O₂ final,
   PHO biofilm mass, φ_b peak, final O₂ and PHO profiles to CSV.
5. **`analysis/PARAMETERS.md`** — full parameter documentation (nominals, sources,
   perturbations, known quirks: baseline sand_pathogen=0, transport_P kills PAT ×40,
   attach_sand hits POM too, transfer rates 60 not 600 — kept for Julia comparability).
6. Smoke test via `matlab -batch`, then launch full campaign in background
   (`analysis/runMatlabCampaign.m`, logs to `SSF/logs/`). → verify: measures.csv rows =
   26, flags OK.

## Deliberate choices

- **Keep `modelLund` transfer rates at 60/d** (not the published 600/300): the MATLAB
  campaign must be comparable to the completed Julia campaign; fixing the preset is a
  separate pending decision.
- Faithful-mirror quirks are preserved, not fixed (they are what Julia ran); documented
  in PARAMETERS.md instead.

## Files touched

- `src/@State/simulate.m` (worktree, matlab-claude)
- `src/presets/pathogenModel.m` (new)
- `analysis/{logOatCampaign,phoBloomStudy,runMatlabCampaign}.m`, `analysis/PARAMETERS.md` (new)
