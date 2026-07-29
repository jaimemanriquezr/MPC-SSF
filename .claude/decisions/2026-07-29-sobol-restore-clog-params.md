# velocity and beta_porosity return to the Sobol design; zeta_0 does not

## The decision

`sobol_sensitivity.jl` sweeps ten parameters. `velocity` and `beta_porosity`,
previously excluded, are restored. `zeta_0`, excluded on the same grounds, stays
out. Both restored parameters use the ranges `morris_screening.jl` uses, and
`beta_porosity` is sampled as β itself rather than as the gap `1 − β`.

## Why

The exclusion comment read: "The clog-driving params (velocity, beta_porosity,
zeta_0) are EXCLUDED — their ranges cause filter failure (undefined QoI)." That
was accurate under the double-well cohesion potential. The corrected potential
removed clogging, and the post-fix Log-OAT confirms it directly: `measures.csv`
in `log_oat_pulse`, `log_oat_startup` and `log_oat_flowstep` reports `flag = OK`
for all 27 parameters, including all three excluded ones.

Leaving the exclusion in place had stopped being merely conservative and become
misleading. The corrected OAT ranks `beta_porosity` **first** in `pulse`
($I_"rms"$ = 0.566) and `velocity` third (0.515). A variance-based study that
omits ranks 1 and 3 of its own top four cannot confirm or refute the co-author
brief's central revision — that biofilm porosity now rivals the attachment
terms — because β is not in the design at all.

`zeta_0` is a different case. It ranks 0.045 in `pulse`, an order of magnitude
below the others, so restoring it would add a design column to obtain a
near-zero index.

Cost is not the deciding factor: `N(k+2)` grows slowly in `k`, so k=10 is 3072
runs against 2560 at k=8, roughly 20% more.

## Verified before committing

A corner test evaluated every parameter at both range extremes with the rest held
at mid-range, on the same proxy the study uses: **0 of 20 evaluations non-OK**,
including `velocity` at 21.6 m/d (3× nominal) and `beta_porosity` at both 0.95
and 0.995. The clogging that motivated the exclusion does not occur.

The corner sweep also shows both restored parameters are live rather than inert:
`velocity` moves mean log-removal from 2.11 to 0.97 across its range, the largest
excursion of any parameter tested.

## β, not the gap

The Log-OAT helper `set_beta_gap()` takes `1 − β`, so reusing it would have been
the smaller diff. It is the wrong choice. Sobol indices are defined against a
distribution over the inputs, and sampling the gap log-uniformly on [0.005, 0.05]
induces a different distribution over β than sampling β uniformly on
[0.95, 0.995]. Reusing the gap helper would have produced indices that are not
comparable with the Morris screen while looking exactly as if they were. A
direct setter costs one line and removes the trap.

## Alternatives rejected

**Keep the eight-parameter design.** Cheapest, and the run was already 14 h in.
Rejected because the result would answer a different question from the OAT
ranking it exists to confirm, and that mismatch would have to be explained in the
paper rather than fixed here.

**Restore all three, k=11.** Rejected on the `zeta_0` reasoning above.

**Run k=8 now and k=10 later.** Two results to reconcile and roughly double the
allocation, to avoid a 20% increase. Rejected.

## Consequence

Jobs `3427791` (k=8, 48 threads, 14.8 h elapsed) and `3428726` (k=8, queued) both
compute the superseded design and were cancelled. The elapsed hours on `3427791`
are lost; it had produced no output in any case, which is what
`.claude/decisions/2026-07-29-sobol-partial-file.md` exists to prevent recurring.

Plan: `.claude/plans/2026-07-29-sobol-restore-clog-params.md`.
