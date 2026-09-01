# MATLAB Log-OAT campaign — parameter documentation (2026-08-18)

Anchor run of the decided methodology in MATLAB (`analysis/logOatCampaign.m`), mirroring
the Julia campaign (`julia/analysis/log_oat_sensitivity.jl`) with the 2026-08-18 design
decisions applied. Plus the PHO bloom study (`analysis/phoBloomStudy.m`).

> **2026-08-27 — everything below the next section documents `Preset = "campaign"`, the
> 2026-08-18 anchor, which ran on the PRE-AUDIT model.** The campaign to be run now is
> `Preset = "workingset"`; its table is the new section
> **"The 29 parameters — `Preset = \"workingset\"`"** immediately after this note. The old
> section is kept unchanged because the completed Julia OAT/Sobol results behind
> `sensitivity.tex` were produced with it and must stay reproducible.

## The 29 parameters — `Preset = "workingset"` (2026-08-27, NOT YET RUN)

Baseline: audited `pathogenModel()` (`.claude/decisions/2026-08-27-pathogen-model-audit.md`)
on the E1–E7 working set, anchored on the **E4 60 d mature filter**
`analysis/probes/data/chain/chain_fld2x_lit_leg6.mat` (`Snapshot=`), pulse scenario,
N = 500, `TPost = 3` d, `MaxDt = 5e-5`, shin scheme, Neumann cohesion BC.
Host configuration: T = 19 °C, ζ₀ = 1, ζ₁ = 0.27, κ = 1e-6, **linear** detachment
0.14·|v_f|/18 ×1, constant liquid transfer ×10, K_DOM 3e-4, K_HPO₄,PHO 1e-6,
η_sand 1500 /m, δ 5 mm, dark floor 0, no respiration reaction.
Influent (kg/m³) = 2× field + the Table B.1 marker in the PAT slot:
`[3.0e-4, 1.0e-3, 0, 5.36e-3, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3]`.
Disturbance: influent PAT ×10 on t ∈ [t_snap + 0.1, t_snap + 0.3) d; measures over
[t_snap + 0.1, t_snap + 3]. Perturbation ×2 / ×½, denominator 2 ln 2, except
`beta_porosity` (constrained scheme, unchanged).

Runner: `slurm/oat_pulse.sbatch` (array 0–28, one parameter per task; each task also runs
the baseline so tasks are independent). Results to `analysis/results/oat/log_oat_pulse_<param>/`.

**Every nominal below is the value the baseline model actually runs at.** That was not true
of the 2026-08-18 table once the kinetics were audited: `logOatCampaign` sets the perturbed
parameter to the *absolute* value 2·nominal or nominal/2, so a stale nominal both perturbs
around the wrong point and silently changes the baseline. Nine nominals moved.

| name | block | nominal (workingset) | was (campaign) | ×2 | ×½ | source |
|---|---|---|---|---|---|---|
| temperature | forcing | 19 °C | 15 °C | 38 | 9.5 | E1–E7 host temperature |
| influent_PAT | forcing | ×1 | ×1 | ×2 | ×0.5 | challenge multiplier on 5.36e-3 (Manriquez B.1) |
| dispersivity | transport | 0.012 m | = | 0.024 | 0.006 | Schijven2013 T1 |
| transport_P | transport | 5.47 /d | = | 10.94 | 2.735 | Lund phase transfer |
| attach_sand | transport | 547 /d | = | 1094 | 273.5 | Diehl2025/Lund b_sand |
| sand_pathogen | transport | 0.06 | = | 0.12 | 0.03 | Schijven2013 T4 α geo-mean (baseline is 0; see quirks) |
| **mu_HET** | kinetics | **2.0 /d** | 1.81e-2 | 4.0 | 1.0 | Reichert2001 k_gro,H,aer |
| **mu_PHO** | kinetics | **2.0 /d** | 5.5 | 4.0 | 1.0 | Reichert2001 k_gro,ALG |
| **d_HET** | kinetics | **0.40 /d** | 2.0 | 0.8 | 0.2 | Wolf2007 b_ina,H |
| **d_PHO** *(new)* | kinetics | **0.276 /d** | — | 0.552 | 0.138 | Campos2006 T3 k_ra |
| **hydrolysis** | kinetics | **3.0 /d** | 0.09 | 6.0 | 1.5 | Reichert2001 k_hyd / Wolf2007 k_h |
| **theta_growth** | kinetics | **dev 0.0725** | 0.047 | 0.145 | 0.036 | Reichert2001 β_H 0.07 → θ 1.0725 |
| **theta_death** | kinetics | **dev 0.08** | 0.066 | 0.16 | 0.04 | Campos2006 θ_kra = 1.08 |
| **K_O2_HET** | kinetics | **2.0e-4** | 3.0e-3 | 4.0e-4 | 1.0e-4 | Reichert2001 K_O2,H = 0.2 g/m³ |
| **K_DOM_HET** | kinetics | **3.0e-4** | 2.0e-4 | 6.0e-4 | 1.5e-4 | **UNCITED** oligotrophic value of the working set; decision pending (`2026-08-26-campos-route-no-respiration.md`) |
| **K_HPO4_HET** | kinetics | **2.0e-5** | 1.4e-8 | 4.0e-5 | 1.0e-5 | Reichert2001 K_HPO4,H = 0.02 gP/m³ |
| **K_O2_PAT** *(new)* | pathogen | **2.0e-4** | — | 4.0e-4 | 1.0e-4 | audited HET row; the manuscript gives no K^O₂_PAT |
| **K_pred** *(new)* | pathogen | 2.0e-3 | — | 4.0e-3 | 1.0e-3 | `thesis_model`; `tab:eco-parameters` gives none |
| marker_growth | pathogen | 0.2 /d | = | 0.4 | 0.1 | Manriquez `results.tex:140` |
| inactivation | pathogen | 0.02 /d | = | 0.04 | 0.01 | Manriquez `results.tex:141` |
| bacterivory | pathogen | 8.0 /d | = | 16.0 | 4.0 | Manriquez `results.tex:142` |
| **zeta_0** | biofilm | **1.0** | 1e2 | 2.0 | 0.5 | E3 working set (`2026-08-25-zeta0-tradeoff.md`) |
| kappa | biofilm | 1e-6 | = | 2e-6 | 5e-7 | interfacial width (sub-grid at this mesh) |
| **zeta_1** | biofilm | **0.27** | 1e-2 | 0.54 | 0.135 | E1 working set |
| detach_scale | biofilm | ×1 | = | ×2 | ×0.5 | scales the **linear** law here, not the campaign √ law |
| light_att_water | light | 0.32 | = | 0.64 | 0.16 | Lund optical depth |
| light_att_sand | light | 1500 | = | 3000 | 750 | Lund optical depth |
| attenuation_P | light | 0.094 | = | 0.188 | 0.047 | Lund self-shading |
| beta_porosity | biofilm | β = 0.99 | = | β = 0.98 | β = 0.95 | constrained scheme, denominator ln 2.5 (decision 2026-08-18) |

Bold = moved from the 2026-08-18 table. *(new)* = added on 2026-08-27.

### Known weaknesses of this design, stated rather than hidden

1. **`K_DOM_HET` has no citation.** 3.0e-4 is the working-set value chosen empirically; the
   Wolf2007 value it replaced (4.0e-3) is an ASM wastewater number. Perturbing an uncited
   nominal produces an uncitable sensitivity.
2. **`marker_growth` may score exactly zero.** With the audited min-Monod set, r7 carries a
   K_HPO₄ term (`pathogen.tex:33`); this campaign runs at influent HPO₄ = 5e-6, so the term
   is limiting but nonzero. At the *manuscript* influent (HPO₄ = 0) it would be identically
   zero. That is why the anchor uses the field influent.
3. **The anchor is the lit arm only.** A dark twin doubles the campaign; the E4 lit/dark
   contrast is 4 % in the bed, so the sensitivity ranking is unlikely to move, but this is
   an assumption, not a measurement.
4. **One 3 d window.** Slow parameters (ζ₀, κ, β) act on the biofilm over weeks; a 3 d
   removal curve cannot resolve them and they will rank low for that reason, not because
   they do not matter.

## Base model — `pathogenModel()` (src/presets/pathogenModel.m)

`modelLund` (physical units, ρ_P = 1117, ρ_L = 998) + Manriquez2026 Table B.4 markers:

| item | value | note |
|---|---|---|
| Lund reactions | µ_HET 1.81e-2, µ_PHO 5.5, d_HET 2.0, d_PHO 0.4, k_h 0.09 /d | inert in flowing phase |
| MarkerGrowth (r7) | 0.2 /d, θ=1.047, Monod O2/NH4/DOM (3e-3, 4e-3, 2e-4) | stoichiometry approximated like HET growth |
| Inactivation (r8) | 0.02 /d, θ=1.08 | |
| Bacterivory (r9) | 8.0 /d, θ=1.08, K_HET = 2e-3 | only reaction active in flowing phase, × water_factor = 1e-3 |
| PAT transport | ×40 the Lund particle rate (5.47 → 218.8 /d) | srun_pathogen convention |
| PAT SandAttachmentFactor | **0.0 at baseline** | perturbed runs use 0.12 / 0.03 — see quirks |
| Cohesion | κ = 1e-6, ζ₀ = **1e2** (overrides Lund 1e6), ζ₁ = 1e-2 | |
| Detachment | 1.4e-5·√(\|v\|/18) | q̂ = 18 m/d hardcoded |
| β (biofilm porosity) | 0.99 | τ (osmosis) = 1e-7 |
| Liquid transfer rates | **60 /d all liquids** | value AT CAMPAIGN TIME, kept for Julia comparability. Fixed to the published 600 (DOM 300) in both presets on 2026-08-18, AFTER the campaign completed — results in this tree predate the fix |

Filter: height 1 m, depth 1 m, sand porosity 0.4, roughness 5e-3, T = 15 °C,
q_in = 7.2 m/d (**fixed — not a sensitivity parameter**), diel light
`0.8·max(sin(2π(t−13/48))+31/50, 0)/(1+31/50)`, attenuation η_water 0.32, η_sand 1500.

Influent (kg/m³): HET 2.68e-3, PHO 1.00e-2, POM 0, PAT 5.36e-3, O2 9.10e-3, IC 6.23e-3,
NH4 2.00e-5, HPO4 0, DOM 1.75e-4 (Table B.1).

## Solver

`simulate` with `TimeStep="adaptive"`, `AdaptiveInitialDt=1e-8`, `AdaptiveMaxDt=3e-6`,
`ImplicitOsmosis=true` (backward-Euler osmosis, new port from Julia — see
implicit-osmosis.md), ncells = 30, CloggingFraction 0.99.

## Scenarios (identical to Julia campaign)

| scenario | initial state | disturbance | window |
|---|---|---|---|
| startup | clean | none (ripening is the disturbance) | [0, 1.5] |
| pulse | mature (3.0 d nominal pre-run) | influent PAT ×10 on t ∈ [0.1, 0.3) | [0.1, 1.5] |
| flowstep | mature (3.0 d) | q_in ×1.4 | [0, 1.5] |

Output: removal L(t) = log₁₀(C_ref/c_PAT,out(t)), 60 frames; measures
I_rms, I_max, D_min, I_min, asymmetry as defined in the Julia driver /
`reports/log-oat-report-2026-08.typ`.

## The 26 parameters

Perturbation is ×2 / ×½ around the nominal (log-symmetric, denominator 2 ln 2) except
where noted. `velocity` is **excluded** (decision 2026-08-18: fixed operating condition).

| name | block | nominal | applied to | source |
|---|---|---|---|---|
| temperature | forcing | 15 °C | filter | Campos2006/Manriquez seasonal |
| influent_PAT | forcing | ×1 | challenge multiplier | Manriquez B.1 |
| dispersivity | transport | 0.012 m | all particles | Schijven2013 T1 |
| transport_P | transport | 5.47 /d | all particles | Lund phase transfer |
| attach_sand | transport | 547 /d | all particles | Diehl2025/Lund b_sand |
| sand_pathogen | transport | 0.06 | PAT SandAttachmentFactor | Schijven2013 T4 α geo-mean |
| mu_HET | kinetics | 1.81e-2 /d | HET growth rate | Campos2006 kgmaxa |
| mu_PHO | kinetics | 5.5 /d | PHO growth rate | Campos2006 kgmaxb |
| d_HET | kinetics | 2.0 /d | HET death rate | Campos2006 T3 kdb |
| hydrolysis | kinetics | 0.09 /d | hydrolysis rate | Campos2006 T3 kh |
| theta_growth | kinetics | dev 0.047 | θ of both growths (1+dev) | Campos2006 T2 |
| theta_death | kinetics | dev 0.066 | θ of both deaths (1+dev) | Campos2006 T2 |
| K_O2_HET | kinetics | 3e-3 | HET growth Monod | Reichert2001 |
| K_DOM_HET | kinetics | 2e-4 | HET growth Monod | Campos2006 T3 ksCd |
| K_HPO4_HET | kinetics | 1.4e-8 | HET growth Monod | Campos2006 T3 ksp |
| marker_growth | pathogen | 0.2 /d | r7 rate | Manriquez B.4 |
| inactivation | pathogen | 0.02 /d | r8 rate | Schijven2013; Manriquez B.4 |
| bacterivory | pathogen | 8.0 /d | r9 rate | Manriquez B.4 |
| beta_porosity | biofilm | β = 0.99 | **constrained: runs at β = 0.98 and 0.95** (gap 0.02/0.05), denominator ln 2.5 | Lund/Melo2005; decision 2026-08-18 |
| zeta_0 | biofilm | 1e2 | cohesion strength | Cahn-Hilliard |
| kappa | biofilm | 1e-6 | interfacial width (sub-grid at this mesh) | Cahn-Hilliard |
| zeta_1 | biofilm | 1e-2 | stable fraction (potential regenerates) | Cahn-Hilliard |
| detach_scale | biofilm | ×1 | detachment prefactor | detachment law |
| light_att_water | light | 0.32 | filter | Lund optical depth |
| light_att_sand | light | 1500 | filter | Lund optical depth |
| attenuation_P | light | 0.094 | all particles | Lund self-shading |

## Faithful-mirror quirks (inherited from the Julia design, preserved deliberately)

1. **Baseline `sand_pathogen` is 0**, not the nominal 0.06 — perturbed runs use ×2/×½ of
   0.06. Both perturbations therefore *increase* removal vs. baseline, which is what
   drives this parameter's asymmetry ≈ 1.7.
2. **`transport_P` overwrites PAT's ×40 transport** — the sweep sets ALL particles to the
   same value, PAT included.
3. **`attach_sand` also hits POM**, whose Lund nominal is 0.
4. β = 0.995 is no longer visited (old gap ×½); the constrained scheme's numbers are NOT
   comparable to the Julia campaign's β row.
5. Liquid transfer rates 60 /d (see base-model table).

## PHO bloom study (`analysis/phoBloomStudy.m`)

Clean start, tsim = 3.0 d (the bloom-development regime), ambient influent, nominal
diel light. Runs: baseline; influent_PHO ×{0.1, 0.3, 3}; mu_PHO ×{0.25, 0.5, 2};
d_PHO ×{0.5, 2, 4}; suppressed_mild (0.3/0.5×/2×); suppressed_strong (0.1/0.25×/4×).
QoIs: max O₂ anywhere over the run ÷ influent O₂ (supersaturation ratio), outflow O₂
(final), min O₂ (final), phototroph biofilm mass ∫(matrix+enclosed)dz, peak φ_b;
per-run final profiles (z, O₂ flowing, O₂ enclosed, PHO biofilm, φ_b) in
`results/pho_bloom/profiles_*.csv`, summary in `summary.csv`.


## Manuscript-experiment suite (`analysis/manuscriptExperiments.m`, 2026-08-19)

Redo of the 10 article experiments + X1 (O2 day/night profiles) on the corrected model
(PG-excess respiration 0.55/d, transfer 600/300). Protocols follow the legacy HEAD
scripts; full table in `.claude/plans/` (manuscript-driver plan). Deviations from the
campaign setup above:

| item | value | why |
|---|---|---|
| kappa | **1e-7** (not the preset 1e-6) | article Table 2 value; preset value is CH-stiff at fine grids |
| detachment, biofilm family (E1–E4) | 0.14·sqrt(v/18) | legacy srun_biofilm.m:13 |
| detachment, pathogen family (E5/E6/E10) | **0.14·sqrt(v/18)** (unified; legacy srun_pathogen.m had 1.4e-5 — with corrected respiration that value clogs the 10-day window) | deviation from legacy, flagged; prefactor explored separately in X2 (`detachmentStudy.m`, sweep 1.4e-5…0.5) |
| respiration Monods | O2 K = 3e-3, NH4 K = 1e-6 | NH4 term is depletion protection (bottom-bed enclosed NH4 crosses zero ~day 5 at 100 cells otherwise); phosphorus removed from respiration entirely (Wolf2007 r6 has none, and K=1e-10 cannot protect explicit steps) |
| light convention | `NormalizedLight=true`: light curve is read as I/I_opt (optimal factor 1.0) | inherited I_opt=1.814e-2 is Wolf's absolute PHOBIA value vs dimensionless curves (peak 0.8) — noon surface sat at ~44x optimal, photoinhibited (decision amendment 2026-08-20b) |
| respiration light response | Monod K/(K+Î), `RespirationLightK=1.0` | replaces 1−Steele complement (wrong shape under normalized light); Wolf's dark switch ≈ K=4.4e-3 |
| grids | 100 cells local, 500 cells cosmos (jobs 3524843/44 MATLAB + 3524845/46 Julia, account lu2026-2-100) | cosmos sbatch needs `#!/bin/bash -l` |
| nominal T | 15 C; seasons 19/3 C | Table 2 (288 K); text's "20 C" not used |

Known HEAD-vs-manuscript-text discrepancies carried (per Jaime 2026-08-19): covered = 1%
(text: 0.1%), scrape depths 0/5/15/25/50 (text: 0/4/8/12), BigPulse = 100x (text: 2x),
inactivation/bacterivory "Experiments 2/3" absent from legacy code (hook kept in E6).
