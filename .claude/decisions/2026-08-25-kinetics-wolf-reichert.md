# Kinetic audit of `modelLund` against Wolf2007 (PHOBIA, Table VI) and Reichert2001 (RWQM1, Table 8) — and the preset fix (2026-08-25)

Requested by Jaime: "Fix the preset. Run an audit on modelLund and Wolf+Reichert for
the growth rates, temperature correction factors and half saturations." Sources read
from `documents/Wolf2007.pdf` and `documents/Reichert2001.pdf` (pdftotext; RWQM1's
exponents survive, PHOBIA's do not — "104" = 1e-4, "6.25" = 6.25e-6 etc., reconstructed
from the cited Henze/Goldman/Eppley values). Model units kg/m³, /d, θ per °C.

## Rates and temperature factors

| process | ours (before today) | Wolf2007 | Reichert2001 | **now** | note |
|---|---|---|---|---|---|
| μ_HET | 1.81e-2 (θ 1.047) | μ_max,H **5.5** (Horn & Hempel 1997) | k_gro,H,aer **2.0**, β_H 0.07 (θ 1.0725) | **2.0, θ 1.0725** | 1.81e-2 was Wolf's I_opt (row slide). Campos2006 0.70–1.0. Reichert chosen: same stack as the algal rows and the O₂ half-sat. |
| μ_PHO | 5.5 (θ 1.047) | none (ETR-driven, q_max,O₂ from ETR_max) | k_gro,ALG **2.0**, β_ALG 0.046 (θ 1.047) | **2.0, θ 1.047** | 5.5 was Wolf's μ_max,H. Campos 1.0–3.0. θ already matched Reichert. |
| d_HET (→POM) | 2.0 (θ 1.066) | b_ina,H **0.4** (Henze) | k_resp,H,aer 0.2 (respiration, O₂-consuming — a different process); no death row for H | **0.4, θ 1.0725** | 2.0 matched Campos k_db avg 2.06 (Jaime 08-20). Wolf's inactivation is the matching process. Campos value stays a calibration alternative. |
| d_PHO (→POM) | 0.4 (θ 1.08) | b_ina,PH **0.09** | k_death,ALG **0.1**, β_ALG (θ 1.047) | **0.1, θ 1.047** | 0.4 was Wolf's *heterotroph* b_ina. |
| respiration_PHO | PG-excess 0.55 | r6: 0.1·q_max, dark-only | k_resp,ALG **0.1**, θ 1.047, K_O₂ 0.2 g/m³ | **0.1, θ 1.047** (decision 2026-08-25-reichert-respiration) | |
| k_hyd (POM→DOM) | 0.09 (θ 1.08) | k_h **3** (Henze) | k_hyd **3.0**, β_hyd 0.07 (θ 1.0725) | **3.0, θ 1.0725** | 33× low; the POM pool was effectively inert. |
| light | Steele, I_opt (normalised 1) | I_opt 1.814e-2 kmol phot m⁻² d⁻¹, K_inh 8e-5 | K_I 500 W/m², exp(1−I/K_I) form | unchanged | |

Wolf2007 has **no temperature dependence** (constant-temperature experiments, pK at 25 °C);
every θ therefore comes from Reichert's β via θ = e^β: β_ALG 0.046 → 1.047, β_H = β_hyd
0.07 → 1.0725. Our old 1.066/1.08 values (Campos) had no Wolf/Reichert basis.

## Half-saturation constants (already corrected this morning; re-checked here)

| constant | ours now | Wolf2007 | Reichert2001 | agreement |
|---|---|---|---|---|
| K_HET,O₂ | 2.0e-4 | 6.25e-6 kmol = 2.0e-4 | 0.2 g/m³ = 2.0e-4 | both ✓ |
| K_HET,DOM | 4.0e-3 (kg COD) | K_S,H,SS 4.0e-3 | K_S,H,aer 2.0 gCOD/m³ = 2.0e-3 | Wolf; Reichert 2× lower |
| K_HET,NH₄ | 1.0e-6 (protective) | 1e-10 kmol ≈ 1.8e-9 | K_N,H 0.2 gN/m³ = 2.0e-4 (2.6e-4 as NH₄) | **sources disagree 10⁵×**; kept Wolf's "N never limits" |
| K_HET,HPO₄ | 2.0e-5 (as P) | — | 0.02 gP/m³ = 2.0e-5 | Reichert ✓ |
| K_PHO,IC | 1.2e-3 (kg C) | 1e-4 kmol CO₂ = 1.2e-3 | — | Wolf ✓ |
| K_PHO,NH₄ | 2.0e-5 | 1.2e-6 kmol = 2.2e-5 | K_N,ALG 0.1 gN/m³ = 1.0e-4 (1.3e-4 as NH₄) | Wolf; **Reichert 6× higher** |
| K_PHO,HPO₄ | 2.0e-5 (as P) | — | 0.02 gP/m³ = 2.0e-5 | Reichert ✓ |
| K_hyd | 0.1 (POM/HET) | K_S,h,X 0.1 | first order in X_S (no half-sat) | Wolf ✓ |
| K_resp,PHO,O₂ | 2.0e-4 | K_S,PH,O₂ 3e-4 kmol = 9.6e-3 | K_O₂,ALG 0.2 g/m³ = 2.0e-4 | Reichert ✓ |

Not present in the model but in both sources: nitrate as an N source (Wolf K_S,PH,NO₃
1.2e-6 kmol; Reichert (9b) with K_N,ALG), phototroph O₂ half-sat on growth, nitrifiers.

## Why Reichert where they conflict

The growth stoichiometry is RWQM1 (9a) to four decimals, the respiration is RWQM1 (10),
and RWQM1 is the only one of the two with temperature factors and a death rate for
algae. Mixing Wolf's rates with Reichert's θ would be a third model. Wolf is kept for
the substrate half-saturations, where RWQM1 either has no value (CO₂) or a coarser one.

## Consequences

- Every run before 2026-08-25 ~21:00 used the old rates; the corrected-K and Reichert-
  respiration arms of today included. Heterotrophs were switched off (μ 0.018 /d) in all
  of them, which is why the O₂ debt of photosynthesis was never paid and HET was 1 % of
  the biofilm. Expect: HET fraction up, POM turnover 33× faster (k_hyd), DOM consumption,
  O₂ supersaturation reduced or gone, phototroph loss 5× slower.
- Manuscript Table 2 must be rewritten with these values and citations.
- SSF.jl `modelLund` must mirror; goldens re-anchored.
- `pathogenModel.m:55` MarkerGrowth (rate 0.2, old K set) still untouched — its own decision.
