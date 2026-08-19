# Phototroph respiration rate: kra = 0.276 /d (Campos2006), constant not dark-only

**Decision:** the new "Phototroph respiration" reaction (plan
`2026-08-18-phototroph-respiration.md`) uses nominal rate **0.276 /d** with
θ = 1.08, runs light-INdependently (constant maintenance respiration), carries the
exact reverse of the phototroph-growth stoichiometry (PHO −1, O₂ −0.9301, IC +0.36,
NH₄ +0.06, HPO₄ +0.01), and is Monod-limited on O₂ (K = 3.0e-3, the Reichert
half-sat already in the model). Enabling it drops the growth reaction's dark floor
(f_dark) from 0.01 to 0.

**Provenance of 0.276 /d:** Campos2006 Table 3, "Algae loss due to respiration and
excretion", kra = 0.0020–0.0210 h⁻¹, average 0.0115 h⁻¹ (source: Brown & Barnwell
1987); ×24 → 0.048–0.504 /d, average 0.2760 /d. θ_kra = 1.08 from the same table
(Temperature correction factors). Extracted 2026-08-18 via pdftotext -layout from
`SSF/documents/Campos2006.pdf`.

**Why constant rather than dark-only:** algae respire in light too; net
photosynthesis = gross − respiration is the standard formulation. A dark-only term
would need a new light-response mode in both solvers (goldens at risk) for less
physiological realism. Constant respiration is preset-only — zero solver changes.

**Why the O₂ Monod on respiration:** shuts the sink off smoothly as water goes
anoxic, preventing negative O₂ without any new guard.

**Alternatives rejected:** r_resp as a fixed fraction of µ_PHO (5% of 5.5 = 0.275/d
— numerically the same, but the literature value carries provenance); fixing the
death stoichiometry instead (death +0.2005 O₂ stays: it is elemental bookkeeping of
PHO→POM, and respiration now carries the metabolic O₂ cost — changing both would
double-count).

**Implemented (2026-08-18), all verified:**
- Julia `modelLund(; phototroph_respiration=0.0)` + forwarding in
  `pathogen_repro.jl` — `Pkg.test` 257 pass, four golden suites untouched.
- MATLAB `modelLund(PhototrophRespiration=0.0)` + `pathogenModel` +
  `phoBloomStudy(Respiration=...)` (matlab-claude worktree).
- Cross-implementation anchor: identical dark run (15 cells, 2e-3 d, implicit
  osmosis) → final O₂/PHO profiles agree to **1.7e-18** (machine precision).
- Dark-column tests in both ports: with respiration on, total O₂ and PHO strictly
  below the floor model's.

**Interesting side observation:** with influent IC = 0, the OLD floor model dies at
the first step ("unphysical concentration in flowing suspension") because dark
floor growth consumes IC from an empty pool — a second concrete demonstration that
f_dark is an artifact term.

**Note for the campaign:** the sensitivity campaign and bloom studies to date ran
with respiration OFF (faithful to the published model). The
`results/pho_bloom_respiration/` study is the first with the split enabled.


## Amendment 2026-08-19 — reformulated to Wolf2007 (PHOBIA) r6, PG-free collapse

Jaime rejected the constant-respiration kinetics after reviewing Wolf et al. 2007
(Biotechnol. Bioeng. 97(5), the PHOBIA model; PDF at SSF/Wolf2007.pdf), reaction r6.
New formulation (his choices: dark-only, no polyglucose pool):

- **Kinetics**: dark-only via Wolf's inhibition switch K_inh/(K_inh + I). New
  `light_inhibition` reaction field in BOTH solvers (Julia `Reaction.light_inhibition`,
  MATLAB `Reaction.LightInhibition`; mutually exclusive with light dependence; 0 = off,
  goldens unchanged). K_inh = 8e-5 kmol(e)/m2/d normalized by I_opt = 1.814e-2 (the Lund
  preset's optimal_light_factor IS PHOBIA's I_opt) -> 4.41e-3 in solver light units.
- **Rate anchor**: 0.1 x q_max (Tillmann & Rick 2001, cited by Wolf2007) = 0.55/d for
  mu_PHO = 5.5 — replaces the Campos kra 0.276/d as the default for studies.
- **Stoichiometry UNCHANGED** (reversed growth row): collapsing r6's polyglucose pool
  into X_PH nets to exactly these coefficients, Y_PH/PG cancels (net biomass -32 kgCOD
  and net N +0.1704 kmol per kmol O2, Y-free; Wolf's net N release = 0.069 kg NH4 per
  kg PHO vs our +0.06 — consistent). Net NH4 UPTAKE is impossible without the PG pool;
  r6's uptake term exists only because PG is N-free carbohydrate.
- **Full r6 (X_PG storage pool + Y_PH/PG ~ 0.63) rejected for now**: new component in
  both ports + golden re-anchoring, for behavior the PG-free collapse already nets out.

Robustness fixes that fell out (MATLAB simulate.m + lookupQuotients.m): no-quotient
models, <2-reaction models, no-light-dependent-reaction models no longer crash;
light normalization guarded when no optimal_light_factor exists (both ports).

Verified: Julia suite 262 pass (incl. new dark-switch tests, goldens untouched);
MATLAB testRespiration.m passes; tiny-model anchor dark AND bright identical across
ports to 10 significant digits. Study: results/pho_bloom_wolf (Respiration=0.55).


## Amendment 2026-08-19 (2) — full r6 with the polyglucose pool

Jaime overruled the PG-free collapse ("It is supposed to consume NH4 and produce PHO —
implement PG"). Both ports now implement the FULL Wolf2007 r6 when
`phototroph_respiration > 0` (plan `2026-08-19-polyglucose-pool.md`):

- New particulate **PG** (CH2O: COD 1.0667, C 0.4, no N/P), 5th particle, transport-
  identical to PHO. Photosynthesis stores f = `pg_fraction` (default 0.2) per unit PHO,
  with O2 +1.0667f / IC -0.4f added to the growth row.
- r6 per unit PG consumed, yield Y = `pg_yield` (default **0.63, ASSUMED** — Wolf2007
  Table VI does not pin Y_PH/PG; ASM heterotroph yield used; candidate for calibration):
  PG -1, **PHO +Y**, O2 -(1.0667-0.9301Y), IC +(0.4-0.36Y), **NH4 -0.06Y**, HPO4 -0.01Y.
  COD and carbon close by construction. Rate unchanged (0.55/d, first-order in PHO,
  dark-only), plus Wolf's quotient Monod PG/(K*PHO+PG) with K_S,PH,PG = 0.005.
- Framework fixes that fell out: MATLAB `lookupQuotients` now handles reactions mixing
  plain and quotient half-saturations (the old code assumed quotient-only dicts).

Verified: Julia 265 tests pass (goldens untouched; PG charges in light, drains in dark);
MATLAB testRespiration.m (r6+PG) passes; cross-port anchor on the charge/discharge run
identical to 10 significant digits (pg_mid, pg_end, o2_end). Study:
`results/pho_bloom_r6full`. The PG-free collapse remains documented above as the netted
form; it is no longer what runs.


## Amendment 2026-08-19 (3) — PG-in-excess variant is the model of record for the O2 goal

Jaime's spec: assume PG is in excess (NOT tracked) and let respiration's light factor be
the exact complement `1 − L_TERM` of the growth Steele factor (activates in dark places —
and, by the same formula, at the photoinhibited surface where growth also idles; flagged
and accepted). Implemented as `modelLund(pg_excess=true)` / `PGExcess=true`, new
`is_light_complement` reaction mode in both solvers (mutually exclusive with the other
light modes; goldens untouched). Stoichiometry = r6 minus the PG column: PHO +Y,
O2 −(1.0667−0.9301Y), IC +(0.4−0.36Y), NH4 −0.06Y, HPO4 −0.01Y — biomass produced, NH4
consumed. Deliberately mass-non-conservative toward the untracked pool. Y stays ASSUMED
at 0.63 (per Jaime); note Y no longer cancels here — the O2/IC coefficients depend on it.

**Outcome across the four variants (baseline effluent O2, influent 9.10 mg/L):**
published 10.47 · r6 no-PG 6.36 · r6 full (tracked PG, f=0.2, Y=0.63) 10.22 ·
**PG-excess 7.44** — goal met: biologically consistent signs AND an effluent deficit;
supersaturation 1.005, bed min 5.31, biomass net-growing (PHO mass 0.357 vs published
0.297). r6-full's excess-return is understood: a tracked pool makes the diel loop
deferred photosynthesis (net O2 = biomass embodied O2); the untracked pool breaks that
identity by design. Anchors: excess variant identical across ports to 10 digits.
Data: results/pho_bloom_{r6full,excess}. Remaining open lever if a still-stronger
bed sink is wanted: heterotroph endogenous respiration (Campos2006 krb avg 1.72/d).
