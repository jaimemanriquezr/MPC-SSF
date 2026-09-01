# Covered/uncovered light insensitivity is a manuscript-level defect, not a tuning gap

## The decision

Escalate the model's near-total insensitivity to light as a **structural defect in the
phototroph -> sand carbon pathway**, to be investigated mechanistically before the
revision. Do NOT close it by tuning light parameters, and do NOT let
`results.tex:184-185` stand as written: it asserts a phenomenon (covered filters need
less scraping) that the model's own figure reproduces only in the third decimal.

Explicitly SEPARATE two claims that the current text conflates:

- **Removal** insensitive to light -- CONSISTENT with the data. Keep.
- **Standing biomass / clogging** insensitive to light -- CONTRADICTED by the data.
  This is the defect.

## Why (measured, with sources)

### What the field data says

`documents/Campos2002.pdf` (Walton WTW, full scale, Bed 9 uncovered / Bed 10 covered,
identical influent, same works -- a controlled light contrast):

| quantity | uncovered (Bed 9) | covered (Bed 10) |
|---|---|---|
| top 0-2 cm sand biomass | max ~60 µg C/g dry sand (>120 at day 56) | -- |
| weighted mean, top 10 cm | logistic, max **22.5 µg C/g at 97 d** (r^2 = 0.98) | **15 µg C/g at 103 d**, ~linear (r^2 = 0.70) |
| depth structure | decreases with depth, significant after 56 d | **no relationship with depth** |
| terminal head loss | reached at ~104 d | **never reached**; still running at 103 d |
| TOC / DOC removal | 25% / 23% | 23% / 23% |

The two decisive rows are the last two. Light produces a **4x difference in standing
biomass and a qualitative difference in clogging**, while leaving removal essentially
identical. Campos attributes the biomass difference to "photosynthetic inputs of
carbon substrates to the sand from the schmutzdecke".

`documents/Demir2017.pdf` corroborates the operational half independently: hydraulic
conductivity falls 13.7 -> 3.28-6.62 um/s over 55 d, monotonically, with "upper layers
decreasing quickly" -- i.e. the operational variable tracks near-surface biomass.

### What the model does

- `manuscripts/AWR-SSF/results.tex:184`: "the uncovered scenario showing slightly more
  growth than the covered one **upon closer inspection**". The Figure 8 bottom panel
  spans phi_b = 0.390-0.405 -- the covered/uncovered difference lives in the **third
  decimal**, against a measured 4x.
- `reports/sensitivity-campaign-2026-08.typ:155-159`: `light_att_sand` ranks 15th-16th
  of 27, `light_att_water` 17th-20th, `attenuation_P` 22nd. "Across x2/x1/2
  perturbations the model output does not respond to them."

So the insensitivity is not a plotting artefact of one figure; it is the model's
measured behaviour across the whole parameter campaign.

### Why this matters more than a 12% discrepancy

The operational claims in the paper are biomass/clogging claims, not removal claims:
`results.tex:166` states the simulation "supports the observation that there is an
increased need to scrape or plough filter beds in summer due to increased clogging from
biofilm formation (Mauclaire2004)". Clogging follows standing biomass -- the quantity
the model gets wrong -- not removal, the quantity it gets right. The paper is therefore
leaning on the model in exactly the place the field data says it is weakest.

## Alternatives rejected

- **Tune the light attenuation coefficients.** Rejected: the sensitivity campaign
  already swept them at x2/x1/2 and the output does not move. There is no coefficient
  value that produces a 4x biomass contrast, because the pathway carrying the effect is
  not attenuation-limited. Tuning would hide the defect and would also invalidate the
  answer already given to the reviewers' light-attenuation concern.
- **Accept it as a model limitation and soften the text.** Rejected for now: Figure 8
  is a covered/uncovered figure and one of the paper's results. Softening the prose
  without understanding the mechanism leaves an unexplained 4x on the table that a
  referee holding Campos2002 will find.
- **Attribute it to the seasonal/temperature work.** Rejected: this is a LIGHT contrast
  at fixed temperature and fixed influent. It is independent of the winter/summer flip
  ([[2026-08-20-winter-summer-flip]]) and must not be folded into it.

## What to investigate

The suspected mechanism is the carbon route Campos names: phototroph growth in the
schmutzdecke -> death -> POM -> hydrolysis -> DOM -> transport into the bed ->
heterotroph growth. Each of those legs exists in the model, so the question is which one
throttles the flux. Candidates, cheapest first:

1. **Light delivered to the schmutzdecke.** With eta_sand = 1500 and (1 - eps_0) = 0.6
   the decay rate is 900 /m -- an e-folding of 1.1 mm -- and light is POINT-SAMPLED at
   cell centres. At N = 100 (dz = 1 cm) the cell-average irradiance in the top bed cell
   is ~9.9x the centre value (sinh(k dz/2)/(k dz/2)). Phototroph production near z = 0
   may simply be under-resolved rather than physically small.
2. **Where the phototrophs are.** Today's runs put 96.9% of biomass in the SUPERNATANT
   ([[2026-08-24-solverA-drift-three-schemes]] follow-up), so the carbon is produced far
   from the sand it is supposed to feed. That is a regression under bisection and could
   be sufficient on its own -- retest this item after it is fixed.
3. **Hydrolysis and transport rates** on the POM -> DOM -> bed leg.

Item 1 is testable without any simulation (it is an arithmetic property of the
attenuation model) and item 2 may dissolve once the regression is found, so BOTH gate
item 3.

## Status

OPEN. No code change made. Recorded so it is not lost behind the numerics work
(regression bisection, zeta_1 recalibration, refinement study), none of which addresses
it.
