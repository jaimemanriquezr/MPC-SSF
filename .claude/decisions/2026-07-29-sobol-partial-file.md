# Partial Sobol indices go to a separate file, not into `sobol_indices.csv`

## The decision

`julia/analysis/sobol_sensitivity.jl` writes incremental results to
`results/sobol/sobol_indices_partial.csv` after each parameter. The canonical
`results/sobol/sobol_indices.csv` is written **only once, on completion**.

## Why

The obvious implementation of "partial writes" is to rewrite `sobol_indices.csv`
after every parameter, so a killed job leaves whatever finished. That is wrong
here, because the file is not scratch space — it holds a **valid N=64 result from
2026-07-22** that is the current best estimate and the one cited downstream.

Rewriting it incrementally means a resubmitted job that dies after 3 of 8
parameters replaces a complete 8-parameter result with a truncated one. The
change intended to prevent data loss would itself cause data loss, and silently:
the file keeps its name, its header, and its plausible shape.

Separating the two files makes the failure mode safe. A run that dies early
leaves `sobol_indices_partial.csv` alongside an untouched `sobol_indices.csv`,
and the difference between "finished" and "partial" is visible from the filename
rather than inferred from a row count.

## Alternatives rejected

**Rewrite `sobol_indices.csv` incrementally.** Rejected above: destroys the
Jul 22 result on any early death.

**Write partials, then rename on completion.** Equivalent safety, but leaves no
partial file after a successful run — so there is no artifact showing the order
parameters were completed in, which is what makes a stalled run diagnosable.
Rejected as strictly less informative for no gain.

**Timestamped output files per run.** Solves the clobbering problem, but every
downstream consumer would need to learn which file is current. Rejected as a
larger change than the problem warrants.

## Evidence

- Pre-existing result: `julia/analysis/results/sobol/sobol_indices.csv`,
  dated 2026-07-22, N=64, 8 parameters.
- Failure this guards against: job `3427791` (N=256, 48 cores, started
  2026-07-29 02:08) ran 14 h with a 66-byte log and no output, because the
  unpatched script writes nothing until the final parameter completes.
- Verified 2026-07-29 with an N=2 run: partial file carried all 8 parameters,
  `y_A.csv`/`y_B.csv`/`y_AB1..8.csv` were all written, and the Jul 22
  `sobol_indices.csv` was restored byte-identical afterwards.

Plan: `.claude/plans/2026-07-29-sobol-partial-writes.md`.
