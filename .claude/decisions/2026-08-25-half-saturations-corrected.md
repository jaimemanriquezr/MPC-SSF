# Half-saturation constants corrected in `modelLund.m` to the audited, source-faithful set (2026-08-25)

## The decision

`src/presets/modelLund.m` (and therefore `pathogenModel`, `probeChain`, every campaign
driver) now carries the half-saturation set of the 2026-08-20 audit:

| reaction | constant | old (published Table 2) | new | source |
|---|---|---|---|---|
| HET growth | K_O₂ | 3.0e-3 | **2.0e-4** | Reichert2001 RWQM1 K_O2,H = 0.2 g/m³ (= Wolf = ASM) |
| HET growth | K_DOM | 2.0e-4 | **4.0e-3** | Wolf2007 K_S,H,SS |
| HET growth | K_NH₄ | 4.0e-3 | **1.0e-6** | Wolf2007 K_S,H,NH₃ = 1e-10 kmol ("N never limits"); 1e-6 is a depletion-protective stand-in |
| HET growth | K_HPO₄ | 1.4e-8 | **2.0e-5** | Reichert2001 K_HPO4,H = 0.02 g P/m³ |
| PHO growth | K_IC | 2.0e-5 | **1.2e-3** | Wolf2007 K_S,PH,CO₂ = 1e-4 kmol |
| PHO growth | K_NH₄ | 1.2e-2 | **2.0e-5** | Wolf2007 K_S,PH,NH₃ = 1.2e-6 kmol (Eppley) |
| PHO growth | K_HPO₄ | 1.68e-4 | **2.0e-5** | Reichert2001 K_HPO4,ALG = 0.02 g P/m³; Campos2006 ksp avg 2.55e-5 |
| Hydrolysis | K_POM/HET | 2.0e-5 | **0.1** | Wolf2007 K_S,h,X = 0.1 kg COD/kg COD |

Requested by Jaime 2026-08-25 ("Correct the modelLund (and ours currently) half
saturations"). Old values retained in a comment in `modelLund.m` for provenance.

## Why

`.claude/decisions/2026-08-20-kinetic-parameter-audit.md` showed every published value
to be an adjacent-row slide of its cited table. The consequence was measured today
(`.claude/CRITIQUE.md` §10): with K_NH₄,PHO = 1.2e-2 the phototroph Monod factor at the
sand surface is 1.7e-4 and growth is 5.4e-4 /d against death at 0.37 /d — phototrophs
cannot grow at any light level. Lit and dark runs were identical at both the manuscript
and the field influent (`chain_z0w_5_n500_{local,dark,field,field_dark}`), the biofilm
was 40 % trapped phototrophs / 58 % POM / 1 % heterotrophs, and K_NH₄,HET = 4.0e-3 gave
heterotrophs a Monod factor of 0.005. The audit's corrected set had existed for five days
only as hand overrides inside individual probes; the preset the campaigns use never
received it.

## Alternatives rejected

- **Campos2006 nitrogen half-sat for phototrophs (ksn avg 0.155 mg N/L → 2.0e-4 as NH₄).**
  Nitrate-based and 10× Wolf; either value restores growth to O(1) (Monod 0.09 vs 0.48
  at the influent). Wolf was kept for internal consistency with the rest of the
  phototroph row, which is PHOBIA throughout. Revisit if calibration wants slower N uptake.
- **Leaving the preset and overriding per probe.** That is what caused the week's results
  to be produced with photosynthesis off. Rejected.
- **Correcting `pathogenModel.m:55` MarkerGrowth (O₂ 3e-3, NH₄ 4e-3, DOM 2e-4 — the old
  HET set) at the same time.** NOT changed: it controls PAT growth and therefore every
  number in the OAT/Sobol campaign behind `sensitivity.tex`. Needs its own decision.
- **`modelLund.m:129` respiration O₂ half-sat 3.0e-3.** NOT changed: it is a switch-off
  term, not a growth limitation; Reichert's 2.0e-4 would be consistent. Flagged.

## Consequences to carry

- **Every result produced before 2026-08-25 evening used the old set**, including the
  deep study, both ζ sweeps, the ζ₀ = 5 recommendation, the 40 d run, the head-loss
  table, and the OAT draft. All must be re-run before use.
- **SSF.jl mirrors `modelLund`** and its golden suites are the cross-implementation
  reference (`SSF.jl/test/golden/`). The Julia preset must receive the same eight values
  and the goldens must be re-anchored, or the two implementations now disagree. Not done
  here (separate repo, Jaime's call).
- `manuscripts/AWR-SSF/results.tex:149-156` (Table 2) must be updated with the new values
  and citations; the audit file gives the source rows.
