# Inflow velocity is a fixed operating condition, not an uncertain parameter

**Decision (Jaime, 2026-08-18):** exclude `velocity` (= `f.inflow_velocity`, q_in,
nominal 7.2 m/d) from the sensitivity analysis. Treat it as a constant of the
operating scenario.

**Why:** the filtration rate is set by the operator, not uncertain in the way
literature-sourced physical/kinetic parameters are. Flow *variability* is already
represented where it belongs — the `flowstep` disturbance scenario (×1.4 surge).

**Consequences:**
- Log-OAT: exclusion is exact post hoc (OAT perturbs one at a time; all other rows were
  computed at nominal velocity). No rerun. Reported design is 26 parameters; velocity's
  raw numbers remain in `julia/analysis/results/log_oat_*/measures.csv`.
  Report updated: `SSF/reports/log-oat-report-2026-08.typ`.
- New leaders: startup — dispersivity 0.745 > sand_pathogen 0.479; pulse/flowstep —
  {beta_porosity, attach_sand, sand_pathogen} ≈ 0.34–0.57.
- Future OAT designs (incl. the MATLAB golden mirror): 26 parameters, 53 runs/scenario.
- Morris/Sobol are NOT exact under post-hoc exclusion (velocity was a design column and
  interacts with attach_sand/sand_pathogen). Their published tables still include it;
  if the fixed-velocity convention is to extend to them, they need a rerun with k−1
  columns. Open question — not decided.

**Alternatives rejected:** keeping velocity as an uncertain parameter (it dominated every
ranking and drowned the physically informative parameters); re-deriving Sobol indices by
conditioning on velocity bins (N=256 too small to condition without new runs).

**Related:** `pathogen_repro.jl:71` hardcodes the detachment reference q̂ = 18 m/d inside
the detachment closure — perturbing velocity never touched q̂. With velocity now fixed,
the q_in/q̂ ratio is constant at 7.2/18. Whether q̂ should equal the nominal operating
rate is a separate unresolved consistency question.
