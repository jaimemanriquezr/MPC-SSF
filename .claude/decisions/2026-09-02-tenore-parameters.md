# Tenore2021: mu_HET = 4.8, nu_P = 52 (k_tot converted to wet mass)

## The decision

From a close reading of `documents/Tenore2021.pdf` (Tenore, Mattei, Frunzo, CNSNS
94:105577, *Modelling the ecology of phototrophic-heterotrophic biofilms*):

1. `mu_HET` 2.00 -> **4.80 /d** (`src/presets/modelLund.m`). Tenore Table 1 `mu_max,2`.
2. `mu_PHO` stays **2.00 /d** -- already Tenore's `mu_max,1`; only the manuscript
   citation was wrong (`Wolf2007` -> `Tenore2021`).
3. `attenuationParticle` (`nu_P`) 0.094 -> **52 m^2/kg**. Tenore's `k_tot` = 210
   multiplied by f_dry = 0.25.

Species mapping confirmed: `f1` = phototrophs, `f2` = heterotrophs (Fig. 1 legend;
`K1,1` is the IC half-saturation for `f1`, `K2,2` the DOC half-saturation for `f2`).

## Why

**The manuscript cited Tenore2021 for a number that is not in it.** `mu_HET` was
1.81e-2 /d. Tenore's heterotroph maximum growth rate is 4.8 /d. 1.81e-2 is that paper's
`Im` = 0.01814 kmol m^-2 d^-1, the OPTIMUM LIGHT INTENSITY -- which `tab:eco-parameters`
already lists separately as `I_opt` = 1.814e-2. One number, two rows, wrong one on the
growth rate.

**nu_P rests on a measurement, not an estimate.** Read off
`analysis/probes/data/chain/chain_fld2x_lit_leg6.mat` (60 d mature anchor), script kept
at `analysis/probes/` conventions, one-off:

| quantity | value |
|---|---|
| particulate concentration, max | 3.61 kg/m^3 |
| particulate concentration, mean in bed | 0.215 kg/m^3 |
| particulate volume fraction, max | 0.0032 |
| AREAL biomass, whole column | **0.221 kg/m^2 wet** = 0.055 dry |
| Tenore biofilm areal (~1 mm at rho = 50) | ~0.05 kg/m^2 dry |

The two systems carry the SAME areal loading. Tenore's coefficient therefore transfers
directly once the wet/dry basis is matched:

| nu_P | eta_P over the column | transmission |
|---|---|---|
| 0.094 (old, Gallegos2000) | 0.021 | **97.9%** |
| 52 (adopted) | 11.5 | 1.0e-5 |
| 210 (Tenore, dry basis) | 46.5 | 6.7e-21 |
| Tenore's own model, for reference | ~10.5 | 4.5e-5 |

At 0.094, biofilm self-shading attenuated light by 2% over the ENTIRE bed -- the
mechanism was effectively switched off. That also reframes the OAT result: `attenuation_P`
ranked near the bottom, which read as "this parameter does not matter" but is better read
as "this parameter was set to a value that disables the mechanism".

Gallegos2000 is not wrong, it is off-context: it measures attenuation by organic matter
suspended in natural water. We apply it to biofilm, which is Tenore's context.

## *** COUPLING TO f_dry -- READ BEFORE TOUCHING EITHER ***

`nu_P = 52` and `f_dry = 0.25` (decision `2026-09-01-f-dry-convention.md`, implement
after Friday) are the SAME conversion applied at different places.

Today the state variable is WET mass (`rho_P` = 1117 is a hydrated-cell density), so
Tenore's dry-basis coefficient is scaled DOWN by f_dry to match. **If f_dry is applied to
`rho_P` and the influents, the state variable becomes dry mass and `nu_P` MUST revert to
210.** Changing one without the other double-counts the water and produces a 4x error in
light attenuation, in whichever direction.

The manuscript row carries the same warning as a comment above it in `results.tex`.

## Alternatives rejected

- **Keep nu_P = 0.094.** Rejected once measured: it makes self-shading a 2% effect,
  which is not a modelling choice anybody made deliberately.
- **Adopt 210 verbatim.** Correct coefficient, wrong basis -- it would multiply wet mass
  by a dry-mass coefficient, over-attenuating ~4x (eta = 46 against Tenore's own 10.5).
- **Adopt Tenore's rho = 50 kg/m^3 as our `rho_P`.** Rejected: different quantity. Cells
  settle in water, so a cell must be denser than water; 1117 vs rho_L = 998 is right for
  a hydrated cell. Tenore's 50 is dry biomass per unit biofilm volume -- a bulk loading
  for a structure that is 90-99% water, consistent with Melo2005 Table 1 (14-91 kg/m^3).
  Both are correct at their own level of description.

## Consequences, NOT yet actioned

Both preset changes invalidate every existing artefact: the E1-E9 chain snapshots, the
OAT campaign E11 (cosmos 3562754 / 3563342 / 3563343), all `rec_*` recreation figures,
and cosmos 3563617 (pat-scraping), which completed today on the old presets. Re-anchoring
is a separate decision and nothing has been re-run.

---

## AMENDMENT 2026-09-02, after cosmos 3564977 arms 0 and 2

**The "self-shading was switched off" reasoning above is wrong and is retracted.** It was
argued from eta_P integrated over the whole 1 m column without checking eta_sand. Measured
decomposition at the final frame of `chain_t22_m20n52_leg3`:

| z | eta_water | eta_sand | eta_P (0.094) | eta_P (52) |
|---|---|---|---|---|
| 0 m | 0.32 | **3.24** | 0.001 | 0.61 |
| 0.01 m | 0.32 | **12.2** | 0.003 | 1.48 |
| 0.05 m | 0.34 | **48.2** | 0.008 | 4.24 |

`nu_sand` = 1500 /m dominates by one to two orders of magnitude at every depth. Light is
already at 2.8% of incident at the top of the bed and 3.5e-6 at 1 cm FROM SAND ALONE, so
the biofilm term cannot matter wherever there is sand.

Measured effect of the change, arm 2 (nu 52) against arm 0 (nu 0.094), both at mu_HET 2.0
under today's code, 30 d lit:

| | areal kg/m^2 | particulate vol.frac. max | HET% | PHO% | PHO @ 1 m |
|---|---|---|---|---|---|
| arm 0, nu 0.094 | 0.2151 | 0.00311 | 65.3 | 33.5 | 10.3 |
| arm 2, nu 52 | 0.2136 | 0.00301 | 65.4 | 33.4 | 10.2 |

**0.7% on areal biomass, 0.1 percentage point on composition.** nu_P = 52 is defensible on
provenance grounds -- Tenore's context is biofilm, Gallegos2000's is suspended matter in
water -- but it is not consequential, and the manuscript should not claim it is.

Two corrections that follow:

1. `attenuation_P` ranking near the bottom of the OAT is CORRECT and well founded, not an
   artefact of a disabled mechanism. The earlier claim in this file that the ranking needed
   reframing is withdrawn.
2. The `ee0e4d3` confound did not materialise. Arm 0 (mu 2.0, nu 0.094, TODAY's code)
   reproduces `chain_fld2x_lit_leg3` (Aug 26) to 4-5 significant figures on every measure
   above, so that 523-line `simulate.m` rewrite was behaviour-preserving for this
   configuration and the old chain is a valid baseline after all.

The mu_HET arms (1 and 3) were still running when this was written; the mu_HET effect is
not yet measured.
