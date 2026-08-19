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
