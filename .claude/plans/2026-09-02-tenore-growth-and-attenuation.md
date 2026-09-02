# Adopt Tenore2021 values: heterotroph growth rate and biofilm light attenuation

## Goal

The manuscript cites Tenore2021 for parameters it does not actually contain. A close
reading (2026-09-02) of `documents/Tenore2021.pdf` (Tenore, Mattei, Frunzo, *Modelling
the ecology of phototrophic-heterotrophic biofilms*, CNSNS 94:105577) gives, in its
Table 1, with `f1` = phototrophs and `f2` = heterotrophs (Fig. 1 legend; `K1,1` is the IC
half-saturation for `f1`, `K2,2` the DOC half-saturation for `f2`):

| Tenore | value | our symbol |
|---|---|---|
| `mu_max,2` (heterotrophs) | 4.8 /d | `mu_HET` |
| `mu_max,1` (phototrophs)  | 2 /d   | `mu_PHO` |
| `k_tot` (light attenuation) | 210 m^2/kg | `nu_P` / `attenuationParticle` |

Done means the presets and the manuscript tables carry these three values, and the
citations on those rows name Tenore2021.

## Steps

1. `src/presets/modelLund.m:82` — `heterotrophGrowth` NominalRate 2.00 -> 4.80.
   Verify: `grep -n "NominalRate=4.80" src/presets/modelLund.m` returns the HET row.
2. `src/presets/modelLund.m:90` — `phototrophGrowth` NominalRate stays 2.00. NO EDIT;
   the working set already equals Tenore's phototroph rate. Verify by inspection.
3. `src/presets/modelLund.m:52` — `attenuationParticle` 0.094 -> 210.
   Verify: `grep -n "attenuationParticle = 210" src/presets/modelLund.m`.
4. `manuscripts/AWR-SSF/results.tex` — `tab:eco-parameters` mu_HET -> 4.80 (cite
   Tenore2021), mu_PHO 2.00 (cite Tenore2021); `tab:general-parameters` nu_P -> 210,
   \unit{\per\kilogram\metre\squared}, cite Tenore2021.
   Verify: `latexmk` still builds, and the three rows read as above.

## Files touched

- `src/presets/modelLund.m`
- `manuscripts/AWR-SSF/results.tex` (on branch `revision`)
- this plan; a decision file follows once the consequences below are settled.

## CONSEQUENCE, not yet accepted by Jaime

Both preset changes invalidate every existing simulation artefact:

- the E1-E9 chain snapshots (`analysis/probes/data/chain/`), which anchor the pulse,
  flowstep and scrape runs;
- the whole OAT campaign E11 (cosmos 3562754 / 3563342 / 3563343, 87 arms);
- every `rec_*` recreation figure in the manuscript;
- cosmos job 3563617 (pat-scraping), which ran on the old presets and completed today.

`nu_P` in particular moves by a factor of ~2234 (0.094 -> 210). It multiplies the areal
biomass density in `eta_P`, so biofilm self-shading becomes ~2000x stronger and light
will not reach below a thin surface layer. Phototroph behaviour throughout the filter
should be expected to change qualitatively, not marginally. `mu_HET` 2.0 -> 4.8 more than
doubles heterotroph growth.

Nothing is re-run under this plan. Re-anchoring is a separate decision.
