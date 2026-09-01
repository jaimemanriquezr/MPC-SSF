# The free-running diagnostic needs its own `matched` branch

## The decision

Add an explicit `matched` branch to the free-running Solver A diagnostic in
`src/@State/simulate.m`, rebuilding Solver B's flux from `uPar`, rather than
letting `matched` fall through to the shin operator.

## Why

The diagnostic block gated the scheme choice on `if ~useBailo`
(`simulate.m:725`, pre-fix). That is a two-scheme test in a three-scheme world:
`matched` is not `bailo`, so it took the `else` path and was advanced by
`lhsPar`, which is assembled from the SHIN stencil
(`CH0 - dt*(S + sparse(D.Rows, D.Columns, D.Values.*lambdaPar))`, line 718).

The comment two lines above the bug already states the requirement — "The
free-running state must use the SAME scheme, or the diagnostic measures the
wrong solver" — so this is the documented invariant being violated, not a new
one being invented.

Consequence if left alone: a three-scheme drift comparison would have reported
`matched` and `shin` as near-identical, and that artefact would have been read
as evidence that matched buys nothing. The wrong conclusion, reached
confidently, from a probe that ran without error.

`vPrev` and `wAdv` are shared with the real state, matching how the shin and
bailo branches already share `S`: the flow field is not part of the CH
discretisation under test — only the operator applied to phi_b is.

Scope: the whole block is inside `if parameters.RecordCflBudget`, so nothing
here reaches `muCH`, `v_b`, or Solver B. Production runs are unaffected.

## Alternatives rejected

- **Run the probe as-is and compare only the one-step metric** (`PhibCH` vs the
  authoritative sum), which *is* scheme-correct for all three because `uCH`
  comes from whichever solve ran. Rejected: the one-step number measures local
  consistency, and the question is whether the mismatch *accumulates*. The 3 d
  data already shows drift peaks creeping upward rather than saturating
  (0.0122 -> 0.0137 -> 0.0149 -> 0.0169), so accumulation is exactly the
  quantity in dispute. Answering with the local metric would dodge it.
- **Drop the matched arm and compare shin vs bailo only.** Rejected: that
  assumes the conclusion. Matched is the only scheme that targets the Solver A/B
  operator mismatch directly, so it is the right control for the claim that
  bailo has absorbed that advantage.
- **Retire `matched` outright without measuring.** Rejected on evidence: the
  numbers previously recorded against it (one-step 1.64e-3 vs shin's 1.08e-4,
  "no drift gain", `simulate.m:645-650`) describe the version with the INVERTED
  sign on `Bop`. The corrected scheme has never been run.

## Evidence

Measured 2026-08-24, cosmos job 3534837 (array 0-2), 20 d, N=100, kappa=1e-6,
MaxDt=3e-6. All three arms Flag=OK and ran an IDENTICAL 6,692,517 steps (dt is
pinned at the MaxDt cap, not CFL-bound, so scheme does not change step count and
wall-time differences are purely per-step cost).

| scheme  | one-step  | accum. max | % of max phi_b | free min   | wall s |
|---------|-----------|------------|----------------|------------|--------|
| shin    | 1.075e-04 | 3.114e-02  | 9.59%          | -7.986e-04 | 4303   |
| matched | 2.220e-16 | 2.615e-12  | 0.00%          |  0.000e+00 | 5086   |
| bailo   | 2.124e-04 | 4.597e-02  | 14.17%         |  0.000e+00 | 6299   |

**The fix took.** matched vs shin free state differs by 9.57e-02 relative; before
the fix it would have been identical to round-off, since both would have been
advanced by `lhsPar`.

**Cross-machine control.** shin run locally (Apple M5) and on cosmos (EPYC 7413)
agrees to 2.1e-13 on `phiA`, 1.3e-13 on `phiS`, over all 6.69e6 steps. Results
are machine-independent; arms from different machines may be compared directly.
Cost ratio local:cosmos = 2371 s : 4303 s = 1 : 1.81.

**Bailo semismooth Newton** (from `rec.bailoNewton`): mean 1.980 it/step,
max 2, non-converged 0%, zero-iter 0.74%. Cost 1.46x shin; matched 1.18x shin.

Raw data: `analysis/probes/data/solverA_drift_<scheme>_n100.mat` (local runs) and
`analysis/probes/data/cosmos/` (the cosmos array, all three).
Figure: `analysis/results/figures/fig8_solverA_drift_n100.png`.
Plan: `.claude/plans/2026-08-24-solverA-drift-three-schemes.md`.

## What this overturned

The hypothesis that motivated the probe -- that bailo's 2026-08-24 v_b fix had
already absorbed matched's operator-consistency advantage -- is WRONG. Bailo has
the WORST accumulated drift of the three (14.17%, against shin's 9.59%). Matching
the exported v_b to the solved mobility did not make bailo's u-row agree with
Solver B's update.

Caveat that must travel with the matched numbers: its ~0 drift is TAUTOLOGICAL.
matched's u-row IS Solver B's discrete update with mu solved implicitly instead
of lagged, so "the difference between Solver A and Solver B" is definitionally
round-off. That is exactly the property component-elimination needs, but it also
means matched provides no independent check on Solver B -- any error in Solver
B's flux becomes invisible rather than showing up as drift.

## Unresolved, and more serious than the drift

The AUTHORITATIVE solutions disagree between schemes:

  shin vs matched   max rel 2.740e-04
  shin vs bailo     max rel 1.557e-02
  matched vs bailo  max rel 1.535e-02

Bailo's phi_b differs from shin/matched by ~1.6% -- a modelling-level difference,
not a diagnostic one, since phi_b is a reported quantity. Which is correct is NOT
established by this probe. Needs a refinement study or an MMS check
(`slurm/mms_ch.sbatch` exists) before any scheme is declared the reference.
