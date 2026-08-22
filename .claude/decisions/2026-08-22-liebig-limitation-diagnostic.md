# Record which substrate binds the Liebig min, as an opt-in solver diagnostic

**Date:** 2026-08-22
**Status:** implemented, verified bit-identical

## The decision

`simulate` gains `RecordLimitation` (logical, default `false`). When set, every frame
additionally records, for the biofilm and flowing regions:

- `Frames.Limitation.Biofilm` / `.Flowing` — `uint8` argmin: which column of `local`
  supplied the min over Monod terms, per cell and reaction (0 = reaction has no Monod terms)
- `.MonodBiofilm` / `.MonodFlowing` — `single`, the value of that min
- `.LightAttenuated` — Î(z,t), the normalised intensity
- `.LightFactor` — the per-reaction light factor actually multiplying each rate
- `.Names` — labels for the argmin values, `[Particles, Liquids, Quotients]`
- `.ReactionNames`

`evaluateReactions` changes signature from `rx = ...` to `[rx, monod, lim] = ...`. MATLAB
permits requesting fewer outputs, so the two remaining single-output call sites are
untouched.

## Why

Two open questions were being answered by *inference from sweep outcomes* rather than by
measurement, and both gate the covered/uncovered driver sweep
(`.claude/plans/2026-08-22-cover-sweep-and-referee-experiments.md`):

1. **Which nutrient limits.** `.claude/decisions/2026-08-20-winter-summer-flip.md` attributes
   the muted light response to phosphorus recycling (influent HPO₄ = 0). But influent
   NH₄ = 2.0e-5 against the published `K_NH4_PHO = 1.2e-2` is a Monod term of 0.0017 — a
   tighter throttle than P. The published and corrected kinetic stacks may therefore be
   limited by *different* nutrients, which decides whether a nutrient add-back sweep relieves
   the right one. Guessing wrong yields a confident null for the wrong reason.

2. **Whether light can reach the bed at all.** `LightAttenuationCoeffSand = 1500` applied to
   `(1−ε₀)(δ/2 + z)` is an effective 900 m⁻¹, implying a euphotic depth of ~5 mm against
   dz ≈ 9.95 mm at 100 cells. Recording Î(z,t) tests that directly instead of by arithmetic.

The deeper reason is structural. Growth is `φ·μ(T)·l·min_j m_j·X`, and a cover of factor s
scales `l` by s *linearly*. Any sub-linear response therefore lives entirely in the feedback,
which is visible only in `min_j m_j`. So the elasticity decomposition

> Δln(growth)/Δln(s) = Δln(l)/Δln(s) + Δln(m\*)/Δln(s)

is measurable, and a total of ≈ 0 with the first term ≈ 1 is the falsifiable signature of the
recycling loop. Without the argmin, that decomposition can only be asserted.

## First results (2026-08-22)

`analysis/probes/probeLimitDiag.m`, 30 cells, 2.0 d, clean start, Table B.1 influent,
manuscript stack (`PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true`):

| reaction | region | phase | binding substrate | mean min-Monod |
|---|---|---|---|---|
| Heterotroph growth | supernatant | Biofilm | **NH4 100 %** | 7.0e-4 |
| Heterotroph growth | supernatant | Flowing | HPO4 70 %, NH4 30 % | 1.7e-3 |
| Heterotroph growth | top 0–2 cm | Biofilm | **NH4 100 %** | 5.8e-2 |
| Heterotroph growth | top 0–2 cm | Flowing | **NH4 100 %** | 1.3e-2 |
| Phototroph growth | supernatant | Biofilm | **NH4 100 %** | 2.4e-4 |
| Phototroph growth | supernatant | Flowing | HPO4 91 %, NH4 9 % | 2.5e-4 |
| Phototroph growth | top 0–2 cm | Biofilm | **NH4 100 %** | 2.0e-2 |
| Phototroph growth | top 0–2 cm | Flowing | NH4 96 %, HPO4 4 % | 4.3e-3 |

**The manuscript stack is nitrogen-limited, not phosphorus-limited**, everywhere except the
flowing supernatant. This corrects the working assumption in
`2026-08-20-winter-summer-flip.md`, which named phosphorus. The 2×2 nutrient add-back
(+P / +N / +NP) in the cover sweep is therefore necessary, not redundant: a P-only sweep
would have relieved the wrong nutrient on this stack.

Same probe, **corrected stack** (μ_HET 1.008, μ_PHO 3.0, source-faithful half-sats), 1.0 d:

| reaction | region | phase | binding substrate | mean min-Monod |
|---|---|---|---|---|
| Heterotroph growth | supernatant | Biofilm | **DOM 100 %** | 4.2e-4 |
| Heterotroph growth | supernatant | Flowing | HPO4 96 %, DOM 4 % | 2.6e-3 |
| Heterotroph growth | top 0–2 cm | Biofilm | **DOM 100 %** | 4.8e-2 |
| Heterotroph growth | top 0–2 cm | Flowing | DOM 87 %, HPO4 13 % | 3.8e-2 |
| Phototroph growth | supernatant | Biofilm | **IC 99 %** | 4.8e-3 |
| Phototroph growth | supernatant | Flowing | **HPO4 100 %** | 6.0e-3 |
| Phototroph growth | top 0–2 cm | Biofilm | NH4 43 %, IC 35 %, HPO4 22 % | **0.703** |
| Phototroph growth | top 0–2 cm | Flowing | HPO4 100 % | 0.311 |

A **third** answer — neither the assumed phosphorus nor the manuscript stack's nitrogen, but
**inorganic carbon** in the supernatant biofilm and **phosphorus** in the flowing phase. This
is the concrete justification for measuring rather than assuming: three plausible priors,
three different answers, all wrong.

### The two blockers are regionally separated

Putting the nutrient and light readouts together for the corrected stack:

| region | nutrient status | light status | why covering does nothing |
|---|---|---|---|
| supernatant | **throttled**, min-Monod ≈ 0.005 (IC / HPO₄) | lit, Î 0.06–0.79 | light is present but growth is nutrient-capped, so removing light changes little |
| top 0–2 cm | **nearly unlimited**, min-Monod 0.70 | Î ≈ 1.4e-27 | nutrients are fine but there is no light to remove |

**Neither blocker alone explains the null; each explains it in a different region.** The
corrected stack largely relieves the bed-top nutrient throttle (0.02 → 0.70 against the
manuscript stack) and thereby *isolates* euphotic depth as the sole remaining bed-side
blocker.

Consequence for the cover sweep: Stage B (nutrient add-back) and Stage B′ (euphotic depth)
are **complementary, not alternatives** — B tests the supernatant limb, B′ the bed limb.
Stage B′ must be promoted from contingent to scheduled. Stage B's add-back levels must also
be re-specified: on the corrected stack the binding species are IC and HPO₄, not NH₄, so a
`+N` arm relieves nothing and should be replaced by `+IC` or by an "all nutrients
saturating" ceiling arm.

Light profile, same run:

```
surface Î max          0.794
Î at z = 0 max         0.0605
Î first bed cell max   1.38e-27
```

**Light is extinguished 27 orders of magnitude within one grid cell of the sand surface.**
Photosynthesis in the bed is not small, it is absent. Campos2002's 0–2 cm sampling layer is
two cells, both of which are optically dark in this model. Any covered/uncovered contrast in
*bed* biomass is therefore mediated entirely by supernatant algae settling and attaching —
it is not a photosynthetic response at all.

This promotes euphotic depth from "fourth candidate blocker" to a measured fact, and it
means Stage B′ of the cover sweep (η_sand ÷ 10, dz ÷ 3) is likely to be needed rather than
contingent.

## Alternatives rejected

- **Infer the limiting nutrient from add-back sweep outcomes.** ~6 core-hours and it confounds
  "the nutrient was not limiting" with "relieving it did not help". The diagnostic costs zero
  compute and is unambiguous.
- **Post-process from saved concentrations.** Would require re-deriving `local`, the
  quotients and the half-saturation lookup outside the solver — duplicated logic that could
  drift from the solver's own. The argmin is already computed inside `min()`; discarding it
  and reconstructing it later is strictly worse.
- **Record unconditionally.** ~13 MB per run and a changed `Results` shape for every existing
  consumer. Default-off keeps every current artefact byte-identical.
- **Log to stdout instead of frames.** Not depth- or time-resolved, and unusable for the
  elasticity decomposition.

## Verification

`analysis/probes/probeLimitDiag.m` plus a three-way bit-exactness check (pre-edit solver via
`git stash`, post-edit with flag off, post-edit with flag on): 30 cells, 0.2 d, all 31
concentration arrays plus the frame times.

```
checked 31 arrays | flags OK OK OK | tFinal 0.200000551814 (all three)
RESULT flag-off bit-identical: 1 | flag-on bit-identical: 1
```

`checkcode` clean. One real bug caught during implementation: the first draft sized the
preallocation from `muRates`, which is built at `simulate.m:239`, *after* the preallocation
block at `:153`. Now uses `numel(model.Reactions)`.

## Files

- `src/@State/simulate.m` — `RecordLimitation` argument; `evaluateReactions` returns
  `[rx, monod, lim]`; new `limitationNames()` local; preallocation, frame write and export.
- `analysis/probes/probeLimitDiag.m` — new, reports argmin occupancy and the light profile
  per reaction × region × phase.
