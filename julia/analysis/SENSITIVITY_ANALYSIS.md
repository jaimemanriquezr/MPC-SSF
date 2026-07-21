# Sensitivity analysis & UQ — setup

Sensitivity analysis of the physical slow-sand pathogen model. Parameter ranges
and the reference SA/calibration methods are drawn from **Campos et al. (2006)**
(biological rates, temperature factors θ, half-saturations) and **Schijven et al.
(2013)** (pathogen attachment / inactivation coefficients). See
`reference_methods.md` for a digest of how each paper did its own SA/calibration.

## Primary method — Logarithmic OAT (`log_oat_sensitivity.jl`)

Per the updated instructions in `Logarithmic_OAT_Sensitivity_SSF_model.pdf`, the
main analysis is a **logarithmic central one-at-a-time** study:

- Build **one mature-filter state** (nominal model run to maturity); every
  perturbed run starts from it, isolating the sensitivity of the *disturbance
  response* rather than of the mature state itself.
- **Disturbance:** an influent pathogen **pulse** (a seeding/challenge event, like
  Schijven's column experiments). Removal `L(t) = log10(C_ref / c_PAT_out(t))`.
- Perturb each θᵢ by **×2 and ×½** (bounded params θ≈1.05, β=0.99 are perturbed on
  their deviation θ−1, 1−β). Measures over the post-disturbance window `[t_d,T]`:
  - `sᵢ(t) = (L⁺ᵢ − L⁻ᵢ)/(2 ln 2) ≈ ∂L/∂ln θᵢ`
  - **`Iᵢ = √⟨sᵢ²⟩`** (RMS log-sensitivity — *primary ranking*)
  - `Iᵢᵐᵃˣ = max|sᵢ|`, `Dᵢᵐⁱⁿ` (effect on worst removal Lmin), `Aᵢ` (asymmetry →
    nonlinearity/threshold).

This is a direct generalization of **Campos (2006)**'s OAT sensitivity coefficient
`SC = COV_out/COV_param` (a normalized log-sensitivity), evaluated on the
pathogen-removal output and the disturbance-response scenario that **Schijven
(2013)**'s scenario-forcing analysis motivates.

Run: `julia --project=julia julia/analysis/log_oat_sensitivity.jl [tmature] [tpost] [ncells]`.
Writes `results/log_oat/{measures.csv, L0.csv, curves_<param>.csv}` and prints
rankings by `Iᵢ` and `Dᵢᵐⁱⁿ`.

## Fuller UQ plan (from `Sensitivity_Analysis_UQ_SSF_model.pdf`)

The log-OAT is stage 1. The broader plan, for when a global/interaction analysis
is warranted:

1. **Numerical uncertainty first** — mesh (Nz = 100/200/400) + tolerance study.
2. **Morris screening** (`morris_screening.jl`, implemented) — global µ\*/σ over the
   same 22 parameters; an alternative/complement to the OAT for interaction effects.
3. **Variance-based Sobol** on the retained 8–12 parameters (`Sᵢ`, `S_Tᵢ`).
4. **Calibration** of an identifiable subset (Bayesian; cf. Schijven's max-likelihood
   + LOO cross-validation, Campos's 300-run Monte-Carlo least-squares).
5. **Uncertainty propagation** to effluent predictions.

Both harnesses take `tsim`/IC as args to screen **start-up**, **mature**, and
**disturbance** scenarios.

## Feasibility — why this runs on a workstation

Each evaluation uses `run_proxy` (`proxy_model.jl`): the **implicit_osmosis**
solver flag keeps τ and β exactly physical but removes osmosis from the CFL,
giving ~20× faster runs with the physics intact (validated to ~1e-4). So a Morris
screen of `r·(k+1)` runs is hours, not cluster-days. The dt is capped at 3e-6
(the faithful ceiling). For the production screen, still prefer a cluster and add
the mesh-convergence check from step 1.

## Quantities of interest (`qois` in `morris_screening.jl`)

| QoI | Meaning |
|-----|---------|
| `Lmean`, `Lmin` | mean / min log₁₀ pathogen removal `log10(c_in/c_out)` |
| `Mb`, `zb` | biofilm mass `∫φ_b dz` and its vertical centroid |
| `pen` | heterotroph (biofilm) penetration depth |
| `o2pen` | oxygen penetration depth |
| `maxHET`, `maxPHO` | peak biofilm heterotroph / phototroph |
| `attPAT` | total attached-pathogen inventory |

## Parameters (22) and range provenance

Ranges are log-uniform (positive) or linear, bracketing the modelLund nominal;
each carries a `source` string citing the paper table + value. Campos rates are
h⁻¹ → ×24 for the model's day⁻¹.

- **Forcing:** `velocity` (Schijven T1), `temperature` (3–25 °C), `influent_PAT`.
- **Transport/retention:** `dispersivity` (Schijven T1), `transport_P`,
  `attach_sand` (Lund b_sand=547), `sand_pathogen` (Schijven T4 sticking α
  5.6e-3–0.71).
- **Biological kinetics (Campos T3/T2):** `mu_HET`, `mu_PHO`, `d_HET` (kdb
  0.05–4.08/d), `hydrolysis` (kh ≈0.067–0.082/d), `theta_growth` (1.047),
  `theta_death` (1.066–1.08), `K_O2_HET`, `K_DOM_HET` (ksCd), `K_HPO4_HET` (ksp
  1e-6–5e-5; **phosphorus is limiting — influent HPO4=0**).
- **Pathogen/marker (Schijven + Manriquez B.4):** `marker_growth` (0.2/d),
  `inactivation` (Schijven µ 0–1.84/d), `bacterivory` (8.0/d).
- **Biofilm mechanics:** `beta_porosity` (0.99), `zeta_0`, `detach_scale`.

Correlations the UQ plan says to preserve at the Sobol stage (not yet enforced):
temperature×irradiation, velocity×dispersivity, growth×half-saturation,
θ×nominal-rate, attachment×detachment, influent nutrient×microbial.

## Extracted reference data

- `data/campos2006/` — 7 CSVs: filtration runs, physical params + **13 θ
  factors** (T2), **34-parameter min/max/avg ranges** (T3), best-10 calibrated
  (T4), headloss calibration (T5). See its `MANIFEST.md`.
- `data/schijven2013/` — 7 CSVs: 13 experiment conditions (T1), MS2 & E.coli
  fitted two-site rate coefficients (T2/T3), schmutzdecke + **sticking
  efficiencies α** (T4), text values, model constants. See its `MANIFEST.md`.

Figure-only data (breakthrough curves, headloss-vs-time, sensitivity coefficients)
was **not digitized** — flagged in each manifest.

## Running

```bash
# smoke: 1 trajectory, short horizon, coarse grid
julia --project=julia julia/analysis/morris_screening.jl 1 0.05 24
# production screen: r trajectories, tsim days, ncells   (cluster recommended)
julia --project=julia julia/analysis/morris_screening.jl 20 0.3 100
```

Writes `results/morris/mu_sigma.csv` (per-QoI µ\*, σ, provenance) and prints a
ranked table. Retain the top-µ\* parameters for the Sobol stage.

## Status

- ✅ **Logarithmic OAT harness** (`log_oat_sensitivity.jl`) — mature-state IC +
  pathogen-pulse disturbance, Iᵢ/Iᵢᵐᵃˣ/Dᵢᵐⁱⁿ/Aᵢ measures, cited nominals. Primary
  method per the updated instructions.
- ✅ Morris screening harness (`morris_screening.jl`) — global µ\*/σ, alternative.
- ✅ Reference SA/calibration methods digested (`reference_methods.md`).
- ✅ Reference data extracted to CSV (Campos2006, Schijven2013); Schijven Fig 1
  measured points digitized to `data/schijven2013/fig1_digitized/`.
- ☐ Mesh-convergence (numerical UQ) study.
- ☐ Sobol variance-based stage on retained parameters.
- ☐ Correlated sampling / Shapley for dependent inputs.
- ☐ Render L⁺ᵢ−L₀ / L⁻ᵢ−L₀ curves (data in `results/log_oat/curves_*.csv`).
