# PAT export non-closure: isolation battery (2026-09-01) — RESOLVED

**Outcome: first-order time-truncation error of the scheme, no defect, no code
change.** Evidence chain and consequences in
`.claude/decisions/2026-09-01-pat-export-closure-resolved.md`. The goal below
was met by the second criterion (irreducible time-splitting error, quantified:
residual halves with MaxDt, O₂ 8.33e-6 → 4.17e-6 → 2.09e-6 over a 4× range).

## Goal

Name the term behind the PAT mass-budget residual of ~5e-5·supply
(2026-08-27 audit, Finding 5) so the log-removal bound can be tightened and the
five blocked pathogen figures drawn. "Done" = one toggle provably kills the
residual (to the 1e-11 level HET/PHO close at), or the residual is shown to be
irreducible time-splitting error with a quantified bound.

## What is already excluded (by reading, this session)

- **Boundary-flux convention mismatch**: `volumeAvgVelocity = q./porosityBoundaries`
  is computed once (simulate.m:136-137), so the solver's outflow flux
  `ε_face(end)·v_avg(end)·g_out` equals the probe's `q·c_out` identically.
- **Implicit dispersion**: `solveImplicitDispersion` (simulate.m:1141) has zero
  flux rows at faces 1, nC, nC+1 and its ε-weighted column sums vanish
  (poroC'·Base = 0), so the tridiagonal solve conserves the ε-weighted stock in
  exact arithmetic. Round-off over 2e5 steps is ~1e-9, two decades below the
  observed residual — but it stays in the battery as an empirical check on the
  reading.

## Finding (2026-09-01, evening): the defect is NOT PAT-specific

Extending the audit to liquids (which the original probe never budgeted) shows
every flowing liquid — O₂, IC, NH₄, DOM — carries the *identical* residual ratio
8.33e-6 of supply (PAT: 5.28e-6; HET/PHO export ~2e-6 of supply, so their clean
residuals were silent, not clean evidence). Attachment, transfer, and every
pathogen-specific term are exonerated; the effect is generic to flowing-phase
transport and proportional to exported mass. Lead hypothesis is now **frame
quadrature of the breakthrough front** (241 frames at 0.0125 d spacing vs a
front crossing the outlet on the same timescale): the budget's trapz of q·c_out
misses a sliver of the front, and the solver itself may be exactly conservative.
Decisive test running: NFrames 241 → 481 → 961; ~4× reduction per doubling
convicts quadrature, a plateau convicts the solver.

## Steps

1. `analysis/probes/probeMassClosureIsolation.m` — probeMassClosure with a
   `Variant` switch and a per-frame residual curve R(t):
   - `baseline` — reproduce the 5e-5 residual → verify: PAT res/supply ~5e-5,
     HET/PHO ~1e-11.
   - `explicitDispersion` — `ImplicitDispersion=false`.
   - `noDispersion` — all particle+liquid Dispersivity = 0.
   - `noTransfer` — PAT TransportRate = 0 (keeps Lund transfer).
   - `explicitOsmosis` — `ImplicitOsmosis=false`.
   - `noOsmosis` — OsmosisRate = 0.
   → verify: exactly one lever (or a named combination) collapses the PAT
   residual; R(t) shows *when* the mass appears/disappears (front passage vs
   steady accumulation).
2. Run all variants locally in background (small N=100 probes, zero biology;
   plots-and-probes stay local per workspace convention).
3. Decision file naming the mechanism + whether a solver fix or a documented
   bound is the resolution — **jointly with Jaime** (this gates the pathogen
   figures and the OAT caveat text).

## Files touched

- `analysis/probes/probeMassClosureIsolation.m` (new)
- `.claude/plans/2026-09-01-pat-export-closure.md` (this file)
- later: `.claude/decisions/2026-09-0X-pat-export-closure.md`
