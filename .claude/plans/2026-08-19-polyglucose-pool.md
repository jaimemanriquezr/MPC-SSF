# Full Wolf2007 r6: polyglucose storage pool (X_PG)

## Goal

Replace the PG-free collapse with the faithful PHOBIA r6 (Jaime, 2026-08-19: "It is
supposed to consume NH4 and produce PHO — implement PG"). Dark respiration becomes dark
GROWTH on internal reserves: biomass produced (+), ammonia consumed (−), polyglucose
consumed, with the pool charged by photosynthesis.

Done = both ports expose the pool behind the existing `phototroph_respiration` opt-in,
goldens untouched (default off), a light→dark test shows PG charging then discharging
with PHO growing and NH4 consumed in darkness, cross-port anchor at golden tolerance,
bloom study rerun (`results/pho_bloom_r6full`).

## Design (mass-unit adaptation of Wolf2007 Tables III/IV)

New particulate **PG** (internally stored polyglucose, CH2O: COD 32/30 = 1.0667 kg O2/kg,
carbon 0.4 kg C/kg, no N/P). Transport-wise identical to PHO (it is inside the cells):
density 1117, dispersivity 1.2e-2, transport 5.47, attachment 547, attenuation 0.094.
Appended as the 5th particle (component order [HET, PHO, POM, PAT, PG | liquids]), so all
Lund-structure index assumptions (POM row 3, reactions 1–5) survive. Influents gain a
PG = 0 slot.

**Photosynthesis stores the f-fraction** (Wolf: storage coupled to growth; f = 0.1–0.4,
default 0.2). Per unit PHO built, additionally +f PG, with the O2/IC of carbohydrate
synthesis added:

| component | coefficient |
|---|---|
| PHO | +1 |
| PG | +f |
| O2 | +(0.9301 + 1.0667·f) |
| IC | −(0.36 + 0.4·f) |
| NH4 | −0.06 |
| HPO4 | −0.01 |

**r6 dark respiration/growth**, per unit PG consumed, yield Y = Y_PH/PG (default 0.63,
the ASM heterotroph yield — Wolf2007's Table VI does not pin it; ASSUMED, revisit):

| component | coefficient |
|---|---|
| PG | −1 |
| PHO | **+Y** |
| O2 | −(1.0667 − 0.9301·Y) |
| IC | +(0.4 − 0.36·Y) |
| NH4 | **−0.06·Y** |
| HPO4 | −0.01·Y |

COD and elemental balance close by construction (each column is the difference of the
PG-synthesis and biomass rows). With Y = 0.63: O2 −0.4807, IC +0.1732, NH4 −0.0378.

**Rate** (Wolf Table IV): 0.55/d (= 0.1·q_max) × PHO × min-Monod{O2 (K = 3e-3),
quotient PG/PHO (K = 0.005 — Wolf's K_S,PH,PG, COD ratio ≈ mass ratio here)} ×
K_inh/(K_inh + I) with K_inh = 4.41e-3 (existing light_inhibition mode). The quotient
Monod uses the framework's existing "PG/PHO" quotient half-saturation support.

## Steps

1. Julia `modelLund.jl`: inside the `phototroph_respiration > 0` branch, add PG,
   f-augmented growth row, r6 row; kwargs `pg_fraction=0.2`, `pg_yield=0.63`.
   → verify: Pkg.test, goldens untouched; light→dark testset.
2. Julia test rewrite: light phase charges PG; dark phase discharges PG, grows PHO,
   consumes NH4 and O2.
3. MATLAB mirror (`modelLund.m`, options PGFraction/PGYield).
   → verify: testRespiration.m rewrite passes.
4. Cross-port anchor: identical light→dark run, profiles at 1e-7 (expect ~1e-15).
5. `phoBloomStudy.m`: influent gains the PG=0 slot when respiration is on.
   → run `results/pho_bloom_r6full` (Respiration=0.55).
6. Amend the respiration decision file; commit both branches + MPCSSF.jl sync + push;
   update the artifact page.

## Files touched

- `julia/src/presets/modelLund.jl`, `julia/test/runtests.jl`
- `src/presets/modelLund.m`, `analysis/testRespiration.m`, `analysis/phoBloomStudy.m` (worktree)
- `.claude/decisions/2026-08-18-phototroph-respiration-rate.md` (amendment)
