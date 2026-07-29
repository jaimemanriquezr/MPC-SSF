# Sobol resume aborts on configuration mismatch rather than falling back

## The decision

`resume=true` reads `results/sobol/run_config.csv` and compares every recorded
field — `N`, `k`, parameter names, `ncells`, `tpost`, `nframes`, `tmature`,
`seed` — against the current run. Any disagreement raises an error and stops.
There is no partial reuse and no automatic fallback to recomputation.

## Why

The cached vectors are QoI evaluations of specific design-matrix rows. `A` and
`B` are drawn from `MersenneTwister(seed)` and their contents depend on `N` and
`k`; the QoI itself depends on `ncells`, `tpost`, `nframes` and `tmature`.

A resume that tolerated a mismatch would combine evaluations of *different*
models into the Saltelli and Jansen estimators. The result is not a visibly
broken run — it is a full set of Sᵢ and S_Tᵢ with plausible magnitudes and
bootstrap intervals, indistinguishable from a valid result once the log scrolls
away. Sensitivity indices are exactly the kind of output nobody re-derives.

Failing loudly costs a resubmission. Failing quietly costs a wrong conclusion
about which parameters matter, which is the entire purpose of the study.

## Alternatives rejected

**Recompute silently on mismatch.** Safe numerically, but it turns a
configuration error into a full-length rerun discovered only from the wall-clock
bill. The user asked for a resume; quietly not resuming is the wrong answer to
give without saying so.

**Reuse only the fields that still match** (e.g. keep `y_A` when just `tpost`
changed). Wrong: `tpost` changes the QoI for every row, so nothing is reusable.
Per-field reuse rules would be a source of exactly the subtle error this guard
exists to prevent.

**Hash the design matrices instead of recording the configuration.** Catches
changes in `A`/`B`, but not in `tpost`/`nframes`/`tmature`, which alter the QoI
without touching the matrices. Strictly weaker.

## Evidence

Verified 2026-07-29 at N=2, `ncells=10`, `tpost=0.1`, `tmature=0.2`:

- Uninterrupted run, then a run resumed after deleting `y_AB4..8` — resulting
  `sobol_indices.csv` **byte-identical** (`cmp`).
- Fully populated directory with `resume=true` — 10 of 10 matrices reused, zero
  evaluations, and zero mature-state builds (the build is lazy, so a complete
  resume does not pay its several-minute cost).
- `ncells` changed 10 → 20 with `resume=true` — aborted with exit 1, naming the
  offending field, and left `sobol_indices.csv` unmodified.

Plan: `.claude/plans/2026-07-29-sobol-resume.md`.
Companion decision: `.claude/decisions/2026-07-29-sobol-partial-file.md`.
