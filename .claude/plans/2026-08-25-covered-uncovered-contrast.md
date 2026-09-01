# Plan: make the covered/uncovered contrast real (uncovered ≫ covered in biomass)

## Goal ("done")

At field influent, ζ₀ = 5, N = 500, 20–30 d, the model reproduces Campos2002's light
contrast qualitatively:

1. **Standing biomass**: uncovered / covered ≥ 2 in the top 2 cm of sand (Campos: ~4×),
   and an uncovered schmutzdecke forms (supernatant mat ≥ 1 mm) while the covered
   filter stays bare.
2. **Removal invariant**: effluent DOC and O₂ removal within ~10 % between the two
   (Campos: 25/23 % vs 23/23 %).
3. **Head loss**: uncovered > covered, both monotone over the run (Campos: uncovered
   reaches terminal head loss ~104 d, covered never); no un-clogging.
4. No effluent O₂ above the influent at any hour (Elemo2024).

Scored by `analysis/results/figures/*_influent_x_light.png`-style comparisons: bed/top-2 cm
integrals, mat thickness, H/H₀(t), effluent O₂ diel amplitude, biofilm composition.

## Where we are (2026-08-25 evening)

| finding | evidence |
|---|---|
| With the published half-sats phototrophs could not grow (Monod 1.7e-4) | CRITIQUE §10; corrected in `modelLund` |
| With corrected half-sats photosynthesis runs (5 g/m³ O₂ diel) but biomass is unchanged: enclosed-phase NH₄ 0.2–1 µg/L vs 16 µg/L flowing (70× gap) — N-transfer/supply limited | `chain_z0w_5_n500_hs_*` |
| O₂ supersaturated because growth/death cycle netted +1.13 O₂ and the PG-excess sink was a biomass source | fixed: RWQM1 respiration (decision 2026-08-25-reichert-respiration) |
| Influent PHO 1e-2 is 15–20× field; at field influent no mat forms lit or dark | `chain_z0w_5_n500_hs_field_*` |
| Model has no nitrate; total inorganic N = 20 µg/L NH₄ (uncited: Chan2018 Table S1 not held) | influent check |

Running now (Reichert respiration): `z0w_5_rr_{lit,dark,field_lit,field_dark,fieldN_lit,fieldN_dark}`
— the 2×2 plus lever L1 below.

## Levers, in the order to pull them (each with the number that motivates it)

**L1 — Nitrogen supply.** Net growth is capped by N: influent N supports ≤ 2.4 g PHO m⁻² d⁻¹.
 a. Stand-in: NH₄_in = 1e-3 (1 mg N/L) at field influent, lit vs dark. RUNNING.
 b. If (a) separates lit from dark: add **NO₃ as a Liquid** and RWQM1 growth (9b):
    NO₃ −0.0600 (as N; −0.266 as NO₃ mass), NH₄ 0, O₂ **+1.2041**
    (8α_C/3 + 8α_H − α_O + 20α_N/7 + 40α_P/31), IC −0.36, HPO₄ −0.01; NH₄-preference
    switch as in RWQM1 (9a/9b split). Both ports + goldens. NO₃_in from Ringsjöverket
    data (ask Jaime; typical lake 0.1–2 mg N/L).
 Verify: enclosed NH₄/NO₃ at the surface no longer < 10 % of flowing at noon; lit bed
 top-2 cm > dark by ≥ 1.5×.

**L2 — Flowing→enclosed transfer of liquids.** `modelLund` liquids `Transport = 600 /d`
(DOM 300). At the surface the biofilm's uptake outruns transfer by ~2 orders (Δc 70×).
 a. Provenance check of 600 /d against Table `rhs-parameters` (`results.tex:98-114`).
 b. Ladder 600 → 6000 → 6e4 at field influent + L1a, lit vs dark, 10 d (cheap: rate
    only enters the source term).
 Verify: enclosed/flowing NH₄ ratio at noon → 0.5+; growth Monod (enc) ≥ 0.3.

**L3 — Phototroph loss rates.** Death 0.4 /d is Wolf's *heterotroph* b_ina (audit);
RWQM1 k_death,ALG = 0.1, Wolf's phototroph 0.09. With respiration 0.1 the total loss
becomes 0.2 /d instead of 0.5. Set death = 0.1 (θ 1.047 per RWQM1). One decision file.
 Verify: lit net growth positive at the surface at noon; dark biomass decays faster
 than lit (it cannot today because both are trapping-dominated).

**L4 — Growth rate.** μ_PHO = 5.5 /d is Wolf's heterotroph μ_max (audit); RWQM1
k_gro,ALG = 2.0 /d, Campos 1–3 /d. Lower it *after* L1–L3, not before: it reduces the
contrast, and its purpose is realism of the O₂ amplitude, not the biomass ratio.

**L5 — Light into the sand.** Campos's 4× is *in the top 2 cm of sand*. η_sand = 1500 /m
(e-folding 1.1 mm) makes the bed dark below one cell; the contrast can then only reach
the sand via carbon export. `coverCases.m` Stage B′ (η_sand 150 /m) already exists.
 Verify: lit − dark difference appears inside z ∈ [0, 2 cm], not only in the mat.

**L6 — Influent** (done): HET 1.5e-4, PHO 5e-4, HPO₄ 5e-6; IC text fixed; POM 0.

Not levers: ζ₀, ζ₁, κ (transport of biomass that has to exist first), attenuation
coefficients at ×2/÷2 (OAT already showed no response).

## Experimental sequence

1. Read the six `rr_*` arms: does Reichert respiration remove supersaturation; does
   L1a separate lit from dark. → decides L1b.
2. L2b ladder (3 pairs × 10 d, local, ~1 h) on top of the best L1 rung.
3. L3 (one pair, 20 d). At this point the criterion 1 should be met or the model has a
   structural gap beyond parameters (then L5, then L1b).
4. Confirm at 30 d with hourly frames; score criteria 1–4; regenerate the manuscript's
   Fig. roofed-results at these settings; write the decision file for the parameter set.
5. Only then: seasons (the same N/loss levers govern the winter/summer ratio) and OAT.

## Files touched

- `src/presets/modelLund.m` (death rate; NO₃ component + (9b) if L1b), `src/presets/pathogenModel.m`
- `analysis/probes/probeChain.m` (options for transfer-rate and death overrides)
- `analysis/probes/coverCases.m` (reuse Stage B′ for L5)
- `.claude/decisions/` one file per adopted lever
- SSF.jl mirror for any preset change (separate repo, Jaime's call)

## Verification

Each rung is a lit/dark pair at identical settings; the score is the ratio table above
plus the enclosed/flowing nutrient ratio at the surface at noon. A rung is accepted only
if lit ≠ dark by more than the run-to-run noise (0.5 % on bed biomass, measured on the
old-K pair) AND removal stays within 10 %. Every accepted rung gets a decision file with
the numbers and the chain-file paths.
