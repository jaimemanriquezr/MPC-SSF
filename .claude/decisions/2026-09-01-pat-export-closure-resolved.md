# PAT export non-closure resolved: first-order time-truncation error, no defect (2026-09-01)

Closes the 2026-08-27 audit's Finding 5 ("the PAT mass budget does not close to
solver tolerance", `2026-08-27-pathogen-model-audit.md`) and the "diagnose
`simulate.m:921-944`" item. Drafted by Claude from the day's runs; Jaime to
confirm the wording before it feeds the manuscript/response letter.

## The decision

**No code change.** The residual is the advection scheme's own O(Δt) truncation
error, visible only for species that actually exit the column. It is documented
and controlled by Δt, not fixed, because there is nothing to fix: the solver is
conservative in the Δt → 0 limit and formally first-order in time.

## The evidence chain (all runs `analysis/probes/probeMassClosureIsolation.m`,
zero biology, N = 100, 3 d, PAT fed; logs `logs/massclosure-iso-*-2026-09-01.log`)

1. **Reproduction**: PAT residual 5.3e-6 of supply vs 6e-11 for HET/PHO —
   but only once PAT is fed (probeMassClosure's PAT influent slot was zero, so
   the original "closed" result never exercised export; HET/PHO export ~2e-6 of
   their supply, so their clean residuals are silent about export-proportional
   error, not evidence against it).
2. **Excluded by argument**: the boundary-flux convention (volumeAvgVelocity is
   static, so the solver's outflow flux is identically q·c_out, simulate.m:136);
   the implicit-dispersion solve (poroC'·Base = 0 — the tridiagonal operator's
   ε-weighted column sums vanish, simulate.m:1141).
3. **Excluded by experiment**: explicit vs implicit dispersion, dispersion off,
   PAT TransportRate = 0 — residual unchanged (~5.3e-6) in all three.
4. **The liquids were the discriminator** (previously never budgeted): O₂, IC,
   NH₄, DOM all show the *identical* residual ratio 8.33e-6 of supply. Nothing
   PAT-specific — attachment, transfer, sand terms all exonerated at once.
5. **Frame quadrature acquitted**: residuals identical to 4 digits at
   241/481/961 frames.
6. **Conviction**: halving MaxDt halves the residual —
   O₂ 8.33e-6 → 4.17e-6 → 2.09e-6 and PAT 5.28e-6 → 2.51e-6 → 1.13e-6 for
   MaxDt 5e-5 → 2.5e-5 → 1.25e-5. Clean first order, converging to zero.

R(t) timing fits: the error is committed while the breakthrough front (the
solution's stiffest transient) crosses the column, then stays put.

## Consequences

- **The P1 log-removal bound stands but is now explained and tunable.** At
  MaxDt = 5e-5 the budget resolution is ~8e-6 of throughput (P1 measured
  8.4e-5 on the 60 d pulse host — same order, stiffer host); each halving of
  Δt buys ~0.3 log of defensible ceiling. The manuscript's 1.48-log claim was
  never at risk; a future >4-log claim would need smaller MaxDt, not a fix.
- **Manuscript/response-letter sentence available**: the marker budget closes
  to first order in Δt (verified over a 4× range), consistent with the scheme's
  formal order; at the operating step the log-removal resolution is ~2.9 log
  (P1 protocol).
- `analysis/testPathogen.m`'s 1e-4 absolute tolerance is justified by
  mechanism, not just by observation — its header can now say why.
- The 2026-08-27 audit's "likely cause: ε·v_avg ≠ q at the outflow face" is
  wrong and superseded (volumeAvgVelocity is static; see 2 above).

## Alternatives rejected

- **Hunting a per-term bug further** (osmosis battery, in-solver per-step
  budget): superseded by the dt-convergence conviction — a genuine leak does
  not converge to zero at O(Δt). The osmosis variants were stopped mid-run;
  noted: explicit osmosis collapses the stable step (>4 h for a run the
  implicit form does in ~10 min), which is why ImplicitOsmosis stays default.
- **Raising the scheme's order / flux-limiting to shrink the constant**: out of
  scope for the revision; PLANS.md territory.
