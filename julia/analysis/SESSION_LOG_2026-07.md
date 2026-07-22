# Session log — pathogen model, fast proxy, and sensitivity analysis (2026-07-21/22)

Reference record of a working session on the Julia port (`MPCSSF.jl`, branch
`julia-port`): filling the physical pathogen model from Manriquez2026, building a
fast evaluation path (implicit osmosis), extracting reference data from Campos2006
& Schijven2013, and setting up + running a sensitivity analysis (Logarithmic OAT +
Sobol + numerical UQ) with a mesh/0-cell investigation.

All work is committed on `julia-port`; commit hashes are given per section. Detailed
companion notes: `SENSITIVITY_ANALYSIS.md`, `reference_methods.md`,
`zero_cell_analysis.md`, `implicit-osmosis` (memory).

---

## 1. Physical marker (pathogen) rates from Manriquez2026 Appendix B  (`38353dd`)

Read Tables B.1/B.4 (rendered as images — `pdftotext` split the value columns) and
filled `pathogen_repro.jl`'s `pathogen_model()`:

| marker reaction | rate (20 °C) | θ |
|---|---|---|
| aerobic growth r7 | 0.20 d⁻¹ | 1.047 |
| inactivation r8 | 0.020 d⁻¹ | 1.080 |
| bacterivory r9 | 8.0 d⁻¹ | 1.080 |

Two corrections: the physical marker has **three** reactions (growth was missing),
and **influent HPO₄ = 0** (Table B.1 — was a guessed 1e-5). Still-uncertain
(flagged inline): K_pred (not tabulated) and the r7 aerobic-growth stoichiometry.

---

## 2. Fast evaluation path — implicit osmosis  (`8633e31`, `e445f23`, `12c1aa5`)

`modelLund` is cluster-stiff (dt~1e-7). Measured cause: the osmosis relaxation
`1/τ = 1e7` is the **sole** binding CFL term — slowing reactions/half-sats does
nothing (51,902 steps either way, dt pinned 9.6e-8).

**`implicit_osmosis` solver flag** (`simulate`, opt-in, default off): the relaxation
`R = (β·φ_B − φ_e)/τ = A − k·φ_W` is **affine in the enclosed-water unknown φ_W**
(k=(1−β)/τ) once φ_matrix, enclosed content, and reactions are lagged explicitly.
Backward-Euler on that stiff linear term → per-cell scalar division (damping
`1/(1+dt·k)`), unconditionally stable, **no Newton** (it is linearly-implicit /
IMEX; the nonlinear couplings are non-stiff and stay explicit). τ and β are left
exactly physical.

Results: 212/212 tests pass with the flag off (clean no-op); **~20× fewer steps**
matching the explicit run to **1.4e-4** at t=0.15. `run_proxy` (in `proxy_model.jl`)
wraps it with `adaptive_max_dt=3e-6` (the faithful ceiling — past it a fast-
depleting species overshoots negative). `modelLund_proxy` adds optional
reaction/half-sat relaxation (lower fidelity).

---

## 3. T = 3 run, frames + profiles  (`6b17cc8`)

`run_proxy(modelLund)` to T=3, 30 cells: 1.2M steps, 69 s, flag OK. Biofilm forms a
**schmutzdecke at the sand surface**, matures to φ_b≈0.16 (modelLund) / 0.33
(pathogen model) by day 3 — no clog. **φ_e = 99·φ_matrix preserved** (osmotic
equilibrium intact under implicit osmosis). 32 frame CSVs + 6 figures in
`results/proxy_T3/`.

---

## 4. Reference data extraction

- **Campos2006** (`data/campos2006/`, 7 CSVs): all 13 temperature factors θ
  (growth 1.047, death-HET 1.066, others 1.08); Table 3 = 34-parameter
  min/max/avg (death kdb 0.0021–0.17/h, hydrolysis kh, half-sats ksCd/ksp/ksn);
  best-10 calibrated sets. Figure-only curves flagged.
- **Schijven2013** (`data/schijven2013/`, 7 CSVs): fitted two-site rate
  coefficients (MS2, E.coli WR1), **sticking efficiency α 5.6e-3–0.71**, katt1
  14–211/d, inactivation 0–1.84/d, 13 experiment conditions.
- **Schijven Fig 1 digitized** (`data/schijven2013/fig1_digitized/`, `7b8a4f5`):
  744 measured breakthrough/seeding points across 8 panels (58 high / 308 medium /
  378 low confidence). Honest limits: tails reliable, day 0–1 peak cloud low-
  confidence, per-experiment attribution best-effort. Fit curves not digitized
  (regenerable from Tables 2–3).
- **SA/calibration methods digest** (`reference_methods.md`): Campos used OAT with
  a sensitivity coefficient SC=COV_out/COV_param + 6-param Monte-Carlo least-squares
  calibration; Schijven used scenario-forcing sweeps + Hydrus-1D/max-likelihood
  with leave-one-out CV.

---

## 5. Sensitivity analysis — Logarithmic OAT  (`d4cee31`, `e517516`, `d378d2f`, `01970f1`, `2de0128`)

Per `Logarithmic_OAT_Sensitivity_SSF_model.pdf` (superseded an initial Morris
harness, `morris_screening.jl`, kept as alternative). `log_oat_sensitivity.jl`:
one shared **mature-filter state** (`final_state` of a nominal proxy run — resume
verified exact), 22 parameters perturbed **×2 and ×½** (bounded θ, β perturbed on
their deviation), disturbance applied, removal `L(t)=log10(C_ref/c_PATout)`.
Measures over the post-disturbance window: `s_i=(L⁺−L⁻)/(2 ln2)=∂L/∂ln θ`, RMS
`I_i` (primary), `I_i^max`, `D_i^min`, asymmetry `A_i`. Runs on the fast proxy.

Two harness fixes during bring-up: 100×→**10× pulse** (a 100× step front is too
sharp for the capped dt → spurious surface clog), and **post-clog frame masking +
clog-driver separation** (a clogged perturbation's zero-filled trailing frames
otherwise inflate I_rms).

### Ranking — pulse scenario (mature + 10× PAT pulse, 30 cells; baseline Lmean=1.34)
| rank | param | I_i | D_i^min |
|---|---|---|---|
| 1 | attach_sand | 0.53 | 0.38 |
| 2 | sand_pathogen | 0.48 | **0.90** |
| 3 | influent_PAT | 0.14 | 0.028 |
| 4 | dispersivity | 0.033 | — |
| 5 | temperature | 0.018 | — |

Clog-driving (separated, log-sensitivity undefined): velocity×2, beta_porosity×½,
zeta_0×½. K_O2/K_HPO4 exactly 0 (P-starved growth — see §8).

### Three-scenario comparison (`log_oat_startup`, `log_oat_flowstep`; `ad7a017`, `4d29e91`)
| param | pulse (mature) | flowstep (×1.4 surge) | startup (clean IC) |
|---|---|---|---|
| velocity | clog-driver | clog-driver | **2.41 (top)** |
| attach_sand | **0.53 (top)** | clog-driver | 0.26 |
| sand_pathogen | 0.48 | **0.34 (top)** | 0.49 |
| dispersivity | 0.033 | 0.008 | 0.39 |
| beta_porosity | clog-driver | clog-driver | 0.23 |
| growth kinetics | ~0 | ~0 | ~0 |

Notes: the ×2 flowstep surge clogged everything (dropped to ×1.4). A parameter's
**role is regime-dependent** (velocity/beta swing between top-sensitivity and clog-
driver; attach_sand leads a mature pulse but destabilizes under flow). The predicted
rise of growth kinetics during startup **did not occur** (P-starvation, §8).

---

## 6. Extensions — Sobol + numerical UQ

- **Sobol** (`sobol_sensitivity.jl`, N=64, 640 runs; `011b41b`): total-effect S_Ti
  **sand_pathogen 0.68 [0.50,0.94] >> attach_sand 0.35 [0.24,0.48] >> all else ~0**.
  The two attachment params own ~all the log-removal variance. Interaction insight:
  sand_pathogen S_Ti−S_i = 0.40 (interacts strongly, ~with attach_sand); attach_sand
  nearly additive (0.05). First-order CIs wide (N=64 first pass; N≥1000 = cluster).
- **Mesh-convergence** (`mesh_convergence.jl`; `3d18f34`): integral QoIs Lmean/Mb
  converge ≲1–2% by 30 cells; Lmin noisier (~5%); **phibmax NOT converged
  (0.25→0.43)**.
- **Mesh-robustness of the sensitivity** (`robustness_mesh.jl`; `72b048b`): top-param
  I_rms ratios 60/30 = sand_pathogen 0.99, attach_sand 1.01, influent_PAT 1.02 →
  the removal ranking is mesh-stable.

---

## 7. The 0-cell and mesh-dependency of clogging  (`b4166e5`, `72b048b`; `zero_cell_analysis.md`)

The mixed 0-cell (centre on z=0, no computational face there) avoids the undefined
flux at the water/sand interface (Diehl2025). Sand attachment is weighted by
`(1−ε)` (`simulate.jl:385`). The regularization scales — roughness **δ=5 mm** and
Cahn–Hilliard width **√κ≈few mm** (κ=1e-6) — are **sub-grid** on any feasible mesh
(dz≈12–50 mm); resolving them needs ncells≳400–2000 (cluster). The schmutzdecke peak
grows ~1/dz (near-singular surface accumulation).

**Rule adopted:** rank on **integral QoIs** (converge + mesh-robust sensitivity),
treat **clogging qualitatively** (peak φ_b and clog times are 1/dz artifacts — which
params drive failure is robust, their timing is not). This is why the clog-drivers
were reported separately and their times demoted to qualitative.

---

## 8. Headline scientific conclusions

1. **Pathogen removal is attachment/retention-controlled, never growth-controlled**,
   confirmed three ways: OAT (attach_sand/sand_pathogen top), Sobol (those two own
   ~100% of variance), and the disproved startup hypothesis (growth kinetics stayed
   ~0 even during ripening). Root cause: **influent HPO₄=0 (Table B.1) starves
   biofilm growth to ≈0 in every regime** — the biofilm builds by filtration of
   influent biomass, not in-situ growth.
2. **Parameter roles are regime-dependent** — a single-scenario ranking would
   mislead; hence the three scenarios.
3. **Numerical validity is bounded by QoI type**: integral removal QoIs are
   converged and mesh-robust (rankings trustworthy); pointwise clog/peak quantities
   are resolution-limited (qualitative only).

---

## 9. Open / deferred (all flagged in the docs)
- Sobol at N≥1000 (cluster) to tighten first-order CIs.
- Correlated sampling / Shapley for dependent inputs.
- Bayesian calibration of the identifiable subset.
- Cluster-scale (ncells≳400) clogging study for quantitative clog timing.
- Section 4 marker r7 stoichiometry + K_pred from Manriquez2026.

## Files (all under `julia/analysis/` unless noted)
Harnesses: `log_oat_sensitivity.jl`, `morris_screening.jl`, `sobol_sensitivity.jl`,
`mesh_convergence.jl`, `robustness_mesh.jl`, `zerocell_diag.jl`,
`render_log_oat.jl`, `render_proxy_T3.jl`. Model: `proxy_model.jl`,
`pathogen_repro.jl`; solver flag in `../src/simulate.jl`. Docs:
`SENSITIVITY_ANALYSIS.md`, `reference_methods.md`, `zero_cell_analysis.md`, this
file. Data: `data/campos2006/`, `data/schijven2013/`. Results: `results/proxy_T3/`,
`results/log_oat*/`, `results/sobol/`, `results/mesh_convergence.csv`.
