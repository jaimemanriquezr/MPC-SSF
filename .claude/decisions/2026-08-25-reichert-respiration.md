# Phototroph respiration re-based on RWQM1 process (10) (2026-08-25)

## The decision

`modelLund` gains `RespirationForm = "reichert" | "pgexcess" | "wolf"` (forwarded by
`pathogenModel`). `"reichert"` implements Reichert et al. 2001 (RWQM1) process (10),
*aerobic endogenous respiration of algae*: first order in PHO, Monod on O₂ with
K_O₂,ALG = 0.2 g/m³ (2.0e-4 kg/m³), light-independent, θ = e^0.046 = 1.047, rate
**k_resp,ALG = 0.1 /d at 20 °C** (RWQM1 Table 6). Stoichiometry is the exact reverse of
the growth row (9a), which our growth row already is:

| | PHO | O₂ | IC | NH₄ | HPO₄ |
|---|---|---|---|---|---|
| growth (9a, unchanged) | +1 | +0.9301 | −0.3600 | −0.0600 | −0.0100 |
| **respiration (10), new** | **−1** | **−0.9301** | **+0.3600** | **+0.0600** | **+0.0100** |

With the RWQM1 algal composition (α_C 0.36, α_H 0.07, α_O 0.50, α_N 0.06, α_P 0.01) the
O₂ coefficient 8α_C/3 + 8α_H − α_O − 12α_N/7 + 40α_P/31 = 0.9300; the 0.9301 in the code
is the same number. `probeChain` now uses `PhototrophRespiration=0.1,
RespirationForm="reichert"`. Requested by Jaime 2026-08-25.

## Why

`.claude/CRITIQUE.md` §10 and the corrected-half-saturation 2×2: with growth un-crippled
the lit filter produced a 5 g/m³ diel O₂ swing (noon effluent 10.1 > 9.1 influent) with no
biomass gain, because a growth→death cycle netted +1.13 O₂ per kg PHO (growth +0.93,
death +0.20) and the only sink — the PG-excess term — is a *biomass source* (PHO +0.63)
whose light switch K/(K+I) at K = 1 leaves it at ~60 % by day. Gross photosynthesis on
recycled N had no matching consumption. RWQM1 closes the loop: respiration consumes the
O₂ growth released, returns the nutrients growth took, and creates nothing from an
untracked pool. The phantom-carbon source (POM export 15.5 > particulate import 12.7 g/m³
at steady state) goes with it.

## Alternatives rejected

- Campos2006 k_ra = 0.276 /d (avg of 0.048–0.504, "respiration and excretion"). Rejected
  in favour of RWQM1 for consistency with the growth row (RWQM1 composition and
  stoichiometry) and with the half-saturation set adopted the same day. Revisit in
  calibration; 0.1–0.28 is the defensible range.
- Fixing the death row instead (make it O₂-consuming). Rejected: RWQM1 death (11) is
  O₂-neutral by construction (Y_ALG,death = 0.62 "set to avoid consumption of nutrients
  and oxygen"); our death row's +0.2005 O₂ is that bookkeeping. Note RWQM1 k_death,ALG =
  0.1 /d against our 0.4 (Wolf's *heterotroph* b_ina, per the audit) — NOT changed here.
- Removing the legacy forms. Kept behind `RespirationForm` so `testRespiration.m`
  (PG pool, dark switch) still exercises them; `""` resolves to the old behaviour.

## Consequences

- Every run before this change used the PG-excess term; the corrected-K 2×2 and the
  ζ₀ = 1 Dirichlet/Neumann pair launched at 19:2x still do (they started before the
  patch). Re-run for any comparison involving O₂ or biomass.
- SSF.jl `modelLund` must mirror `RespirationForm` and the goldens re-anchored.
- Manuscript: the "dark respiration term" explanation Jaime owes the reviewers now
  describes RWQM1 (10), not Wolf r6 / PG-excess.
