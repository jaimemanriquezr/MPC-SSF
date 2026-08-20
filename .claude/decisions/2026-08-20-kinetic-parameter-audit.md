# Kinetic-parameter audit: growth rates and half-saturation constants (2026-08-20)

## The finding

Manuscript Table 2 (results.tex:135-155) / modelLund.m kinetic constants were audited
against their own citations (Wolf2007 Table VI, Campos2006 Table 3, Reichert2001 RWQM1
— all in SSF/documents/). Every growth rate and every half-saturation constant is wrong
against its cited source. The recurring failure mode is a "row slide": adjacent source-
table rows shifted by one position, mantissas surviving with mangled units/exponents.

## Growth rates

| ours | value | cited | source actually says |
|---|---|---|---|
| mu_HET | 1.81e-2 /d | Tenore2021 | numerically Wolf's I_opt = 1.814e-2 (adjacent row); Wolf mu_max,H = 5.5/d, Tenore-group photogranule model 4.8/d, Campos bacteria 0.70-1.0/d, ASM 6/d |
| mu_PHO | 5.50 /d | Wolf2007 | exactly Wolf's mu_max,H (heterotrophs, Horn & Hempel); Wolf phototrophs have NO mu_max (ETR-driven); Campos algae 1.0-3.0/d, Reichert k_gro,ALG = 2.0/d |

Jaime's correction (2026-08-20): mu_HET = 0.042*24 = 1.008/d, mu_PHO = 0.125*24 = 3.0/d
(Campos2006 Table 3 upper range).

## Half-saturation constants (model units kg/m3)

| ours | value | cited | source-faithful value | note |
|---|---|---|---|---|
| K_HET_O2 | 3.00e-3 | Reichert2001 | 2.0e-4 (RWQM1 K_O2,H,aer = 0.2 gO/m3 = Wolf = ASM) | 15x high |
| K_HET_NH4 | 4.00e-3 | Wolf2007 | ~1.7e-9 (Wolf K_S,H,NH3 = 1e-10 kmol, "N never limits") | ours = Wolf's K_S,H,SS substrate row; THE heterotroph crippler (Monod 0.005 at influent) |
| K_HET_HPO4 | 1.40e-8 | Reichert2001 | 2.0e-5 (K_HPO4,H,aer = 0.02 gP/m3) | 1400x low, benign |
| K_HET_DOM | 2.00e-4 | Wolf2007 | 4.0e-3 (K_S,H,SS) | ours = Wolf's K_S,H,O2 converted; (NH4,DOM) pair = Wolf's (SS,O2) slid |
| K_PHO_IC | 2.00e-5 | Wolf2007 | 1.2e-3 kgC/m3 (K_S,PH,CO2 = 1e-4 kmol) | ours = Wolf's NH3 half-sat correctly converted, landed one row off |
| K_PHO_NH4 | 1.20e-2 | Wolf2007 | 2.0e-5 (1.2e-6 kmol NH3) | mantissa kept, conversion botched 600x; THE phototroph crippler (Monod 0.0017) |
| K_PHO_HPO4 | 1.68e-4 | Reichert2001 | 2.0e-5 (K_HPO4,ALG = 0.02 gP/m3; Campos ksp avg 2.55e-5 agrees) | 8x high |
| K_hyd | 2.00e-5 | Wolf2007 | 0.1 kgCOD/kgCOD (K_S,h,X) | 5000x low -> hydrolysis unsaturated |

## Why it matters

The two "cripplers" explain why the 2026-08-20 growth-rate probe (probeRates.m) changed
nothing: heterotroph and phototroph growth are half-sat-limited, not rate-limited. With
source-faithful values the heterotroph DOM->O2 sink (stoich -1.2317) opens (~200x release)
and in-filter photosynthesis un-cripples (~300x) — the O2 supersaturation anomaly that
the phototroph-respiration feature has been compensating for is plausibly rooted here.

## Corrected probe set (probeHalfSats.m, analysis/probes/)

HET: O2 2.0e-4, DOM 4.0e-3, NH4 1.0e-6 (protective stand-in for Wolf's ~0), HPO4 2.0e-5.
PHO: IC 1.2e-3, NH4 2.0e-5, HPO4 2.0e-5. Hydrolysis POM/HET 0.1.
Plus corrected growth rates. Results in probe outputs (halfsat_{noresp,resp}.mat).

Status: investigation only — model of record unchanged pending co-author decision.

## Addendum: death/inactivation rates (2026-08-20, evening)

Prompted by Jaime's question on inherent algal death in the sources:
- Campos2006 has NO separate algal death rate: k_ra (respiration/excretion) IS the
  entire inherent algae loss (Eq. 1 supernatant, Eq. 23 bed) plus protozoa grazing.
  Bacteria do get k_db = 0.05-4.1 /d, avg 0.086/h = 2.06 /d.
- Wolf2007 Table VI: inactivation b_ina,PH = 0.09 /d (assumed), b_ina,H = 0.4 /d
  (Henze). Inactivation -> inert biomass, no O2 release.
- Our d_PHO = 0.4 (cited Wolf) = Wolf's HETEROTROPH b_ina,H — fourth row slide,
  same phototroph-takes-heterotroph pattern as mu_PHO. Wolf's phototroph value: 0.09.
- Our d_HET = 2.0 (cited Wolf) matches Campos's avg k_db = 2.06 — right magnitude,
  wrong citation.
- Consequence for the endogenous proposal: the double-counting caveat resolves.
  Source-faithful loss structure = death 0.09 (PHO) / ~2.0 (HET) with O2-NEUTRAL
  stoichiometry (biomass -> POM/inerts; our current death rows release O2, matching
  neither source) + endogenous respiration 0.276/1.72 carrying all O2 demand.
