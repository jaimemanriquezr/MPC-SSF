# Phototroph loss = Campos k_ra; no separate respiration reaction (2026-08-26)

## The decision

`modelLund` "Phototroph death" now has NominalRate **0.276 /d**, θ **1.08** — Campos2006
Table 3 k_ra, "algae loss due to the combined effects of respiration and excretion"
(Brown & Barnwell 1987; range 0.048–0.504 /d). It is the single phototroph loss process,
as in Campos2006 (Eqs. 1, 23) and Wolf2007 (b_ina,PH). The RWQM1 respiration reaction
added on 2026-08-25 is no longer instantiated by default: `probeChain` runs
`Respiration = 0` and zeroes the growth reaction's dark floor itself; `RespirationForm`
stays in the preset so the 2026-08-25 arms remain reproducible. Jaime's choice of 0.276
over Reichert's death + respiration = 0.2 (2026-08-26).

## Why

- The respiration sweep on the final preset (`chain_z0w_5_fp_lit_r{0,0p276,0p5,1}`)
  showed effluent O₂ is insensitive to the phototroph respiration rate once
  heterotrophs are active (10 d: 3.9 / 1.6 / — / 1.05 g/m³ for r = 0 / 0.1 / 1.0; no
  value supersaturates, none reaches the observed 3–8 window at manuscript loading;
  the bed's demand on the particulate load dominates). Respiration only sets the
  standing phototroph mass.
- One loss term is what both cited SSF/biofilm models do and removes the "new dark
  respiration term" the reviewers would otherwise need explained.

## What it costs (the O₂ accounting, to be stated in the manuscript)

Death keeps RWQM1 (11) stoichiometry: PHO −1, POM +0.6316, O₂ +0.2005, NH₄ +0.0221,
HPO₄ +0.0037. Per kg of PHO that dies the O₂ debt is paid downstream: 0.63 POM →
0.63 DOM (k_hyd 3 /d) → 0.40 HET consuming 0.49 O₂; net −0.29 vs +0.93 released on
growth. A growth → death cycle therefore nets +0.64 O₂ (RWQM1 respiration made it 0).
Supersaturation is possible again wherever a lit mat grows on recycled nutrients and
heterotrophs cannot consume the resulting DOM — which is the case today at field DOM:
K_DOM = 4 mg/L (Wolf K_S,H,SS, an ASM wastewater value) gives μ_eff = death at 1 mg/L,
so 20 d at DOM_in = 1 mg/L removed 2.7 % (`chain_z0w_5_fp_fieldDOM_*`). **This decision
makes an oligotrophic K_DOM (0.05–0.5 mg/L, AOC-type kinetics) obligatory, not optional.**
Separate decision pending.

## Alternatives rejected

- Reichert death + respiration merged at 0.2 /d, θ 1.047: consistent with the rest of
  the algal row, but Jaime chose the SSF-specific Campos value.
- Keeping RWQM1 (10): closes the O₂ cycle exactly but adds a process neither SSF model
  has and that the sweep showed to be inert for the effluent.

## Consequences

- All arms before 2026-08-26 ~00:30 used either PG-excess or RWQM1 respiration.
- SSF.jl mirror + goldens; manuscript Table 2 row "phototroph death" → k_ra with
  Campos citation; the respiration paragraph drafted for the reviewers is withdrawn.
