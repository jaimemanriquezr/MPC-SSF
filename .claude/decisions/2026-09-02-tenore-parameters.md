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
| phi_b, max | 0.0032 |
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
