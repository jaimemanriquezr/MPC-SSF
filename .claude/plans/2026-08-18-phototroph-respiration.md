# Phototroph respiration: split growth into photosynthesis (O₂+) and respiration (O₂−)

## Goal

Give phototrophs a real oxygen cost. Today their entire metabolism is ONE reaction —
"Phototroph growth" (rate 5.5/d, θ=1.047, Steele light response `I = max(f_dark,
I_eff·e^(1−I_eff))`, stoichiometry PHO +1, O₂ +0.9301, IC −0.36, NH₄ −0.06, HPO₄ −0.01)
— so in darkness the floor f_dark = 0.01 makes algae weak O₂ *producers*, and with
f_dark = 0 they are inert. Real dark algae are a *sink*. The feature splits the
metabolism into:

1. **Photosynthesis** — the existing reaction, light-dependent, `f_dark = 0` (the floor
   retires; gross photosynthesis is genuinely zero in darkness).
2. **Respiration** — a NEW first-order reaction in PHO, light-INdependent (runs day and
   night, as maintenance respiration does), with the growth stoichiometry reversed:
   PHO −1, O₂ −0.9301, IC +0.3600, NH₄ +0.0600, HPO₄ +0.0100.

Reversing the photosynthesis stoichiometry makes the pair elementally consistent by
construction: respiration undoes exactly what photosynthesis built, so no new elemental
bookkeeping is introduced. Net behavior: in light, net O₂ = 0.9301·(µI − r_resp)·PHO;
in darkness, −0.9301·r_resp·PHO — the correct sign in both regimes.

Done = both ports expose the split behind an opt-in flag, goldens untouched, a
dark-column test shows O₂ *decreasing*, and the PHO bloom study rerun under the split
shows bed O₂ drawdown and outflow O₂ ≤ influent for the nominal bloom.

## Evidence this is the right target (2026-08-18 studies, MATLAB anchor)

- Bloom study (f_dark = 0.01): outflow O₂ 10.5 vs influent 9.1 mg/L; d_PHO ×4 removes
  55% of algae yet RAISES O₂ to 11.3 — wrong sign, no respiration cost anywhere.
- f_dark = 0 rerun (`results/pho_bloom_dark0/`, complete, 12/12 OK): the dark floor
  explains only ~10% of the anomaly — baseline outflow drops 10.466 → 10.324 mg/L
  (excess over influent 1.366 → 1.224); the d_PHO wrong sign persists in full
  (9.89 → 11.17 mg/L across the d sweep with the floor off). The supernatant-produced
  excess advects through the bed because nothing consumes it.
  **⇒ only a genuine dark sink changes the bed profile — i.e., this feature.**

## Design

### New reaction (both ports, name "Phototroph respiration")

| field | value | rationale |
|---|---|---|
| order | PHO → 1.0 | first-order in phototroph biomass |
| nominal_rate `r_resp` | **0.05–0.10 × µ_PHO** as default (0.275–0.55 /d); extract Campos2006 kr (Table 2/3) during implementation — our θ_death source note already cites "θ_kr* = 1.08", so Campos carries a respiration constant | typical algal maintenance fraction; literature-anchored |
| temperature_correction_factor | 1.08 (Campos θ_kr) | |
| is_light_dependent | **false** | maintenance respiration runs in light too; net photosynthesis emerges as gross − respiration. Avoids any solver change — no "inverse light response" mode needed |
| half_saturation_constants | O₂ → 3.0e-3 | Monod on O₂ so respiration shuts off as the water goes anoxic — prevents negative O₂; reuse the Reichert O₂ half-sat already in the model |
| stoichiometric_coefficients | PHO −1, O₂ −0.9301, IC +0.36, NH₄ +0.06, HPO₄ +0.01 | exact reverse of photosynthesis |
| efficiency_flowing | 0.0 in the pathogen model (matches every other Lund reaction) | |

### Photosynthesis reaction changes

- `minimum_light_factor`: 0.01 → **0.0** (the floor's job is taken over by respiration).
- Everything else unchanged.

### Opt-in switch — goldens must not move

- Julia: `modelLund(; phototroph_respiration=0.0)` — 0.0 means "reaction absent"
  (or not appended), so the four golden suites and every existing script reproduce
  bit-identically. `pathogen_model()` and `modelPathogen()` forward the kwarg.
- MATLAB: `modelLund(PhototrophRespiration=0.0)` and `pathogenModel(...)` mirror it.
- The oxygen/bloom studies and future campaigns pass the literature value explicitly.

## Steps

1. **Extract Campos2006 kr** (documents/Campos2006.pdf, Tables 2–3) → decide the default
   `r_resp` and record it with provenance in the presets' comments.
   → verify: value + unit conversion (×24 h→d) written in a decision file.
2. **Julia**: add the reaction to `modelLund.jl` behind the kwarg; forward from
   `modelPathogen.jl` and `pathogen_repro.jl:pathogen_model`.
   → verify: `Pkg.test` — all four golden suites still pass with the default off.
3. **MATLAB** (matlab-claude): same in `modelLund.m` + `pathogenModel.m`.
   → verify: pulse-scenario 3-param smoke reproduces the current measures.csv rows
   exactly with respiration off.
4. **Unit test (both ports)**: dark column (light_irradiation ≡ 0), seeded PHO, no
   influent O₂ → O₂ strictly decreases and floors at ~0 (Monod), never negative;
   IC rises by 0.36/0.9301 of the O₂ consumed.
   → verify: new test in `julia/test/runtests.jl` + a MATLAB assert script.
5. **Cross-implementation anchor**: one nominal bloom-study run with identical
   respiration settings in both ports; compare final O₂/PHO profiles at atol 1e-7
   (the pathogen-golden tolerance).
   → verify: comparison script output committed alongside.
6. **Science rerun**: `phoBloomStudy(DarkRespiration=0, Respiration=r_resp)` →
   `results/pho_bloom_respiration/`. Expect: O₂ declining through the bed, outflow ≤
   influent at nominal bloom, and d_PHO regaining the physical sign (more death → less
   algae → less total O₂ turnover).
   → verify: summary.csv sign checks; add the third line to the oxygen-mechanism
   hierarchy (floor=0.01 / floor=0 / respiration) in the reports and the artifact page.

## Alternatives rejected

- **Dark-only respiration** (rate × (1 − light response)): physiologically weaker (algae
  respire in light too) and requires a new light-response mode in BOTH solvers — a solver
  change with golden implications, for less realism. Constant maintenance respiration is
  the standard formulation (ASM/River-Water-Quality-Model lineage) and is preset-only.
- **Fixing the death stoichiometry instead** (death −O₂): death's +0.2005 O₂ is elemental
  bookkeeping of PHO→POM conversion, defensible as written once respiration carries the
  metabolic O₂ cost. Changing it would double-count. Leave death alone; note in the
  manuscript that decay-to-POM is O₂-neutral-ish by design.
- **Reusing f_dark as the respiration knob**: it multiplies the *photosynthetic*
  stoichiometry, so it can never have the right sign. Retire it (also fixes the
  misleading `results.tex:157` table row "growth-rate factor in dark respiration").

## Files touched

- `julia/src/presets/{modelLund,modelPathogen}.jl`, `julia/analysis/pathogen_repro.jl`
- `src/presets/{modelLund,pathogenModel}.m` (matlab-claude)
- `julia/test/runtests.jl` (+ small MATLAB test script)
- `analysis/phoBloomStudy.m` (new Respiration option)
- Manuscript later: `ecomodel.tex` reaction table, `results.tex` parameter table (replace
  f_dark row), oxygen discussion in `results.tex` — separate task on Overleaf.

## Open questions for Jaime

1. Default `r_resp`: take Campos2006 kr as-is, or as a fraction of µ_PHO? (Plan assumes
   Campos wins if the two disagree.)
2. Should heterotrophs get the same treatment later (their death is also O₂-positive,
   +0.0234)? Out of scope here, but the same pattern applies.
3. Does the split go into the Manriquez2026 revision itself, or into the response as
   "identified fix, future work"? (Deadline Sept 1 argues for the latter unless the
   rerun is quick — it is: the bloom study is 22 min.)
