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
