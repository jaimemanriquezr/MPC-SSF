# Sobol: flushed progress and partial writes

## Goal

`julia/analysis/sobol_sensitivity.jl` must never lose a whole run to a wall-clock
kill. Two defects, both live in job `3427791`:

1. **Silent log.** Julia buffers stdout when Slurm redirects it to a file. After
   14 h the log was 66 bytes, so there was no way to tell "working" from "stuck"
   without `sstat`.
2. **All-or-nothing output.** `sobol_indices.csv` is written once, after the last
   parameter. A kill at any earlier point loses every completed evaluation.

Done means: a job killed at the wall leaves behind the QoI vectors for every
completed design matrix and the indices for every completed parameter.

## Steps

1. **Flush after every progress print.**
   Verify: run locally with stdout redirected to a file; `wc -c` on that file is
   non-zero before the run finishes.

2. **Persist each design-matrix QoI vector as it completes** — `y_A.csv`,
   `y_B.csv`, `y_AB<i>.csv` in `results/sobol/`.
   Verify: after a small local run, those files exist and `y_A.csv` has N rows.

3. **Write indices incrementally to `sobol_indices_partial.csv`** after each
   parameter; write the final `sobol_indices.csv` only on completion.
   Verify: interrupt a local run mid-loop; the partial file holds the finished
   parameters and the pre-existing `sobol_indices.csv` is untouched.

4. **Print elapsed minutes** on each milestone so progress is readable from the
   log alone.
   Verify: log lines carry increasing minute stamps.

## Why the partial file is separate

`results/sobol/sobol_indices.csv` currently holds a **valid N=64 result from
Jul 22**. Writing partial rows straight into it would destroy that result the
moment a resubmitted job died early — trading one failure mode for a worse one.
The real output file is replaced only by a run that finished.

## Files touched

- `julia/analysis/sobol_sensitivity.jl`

## Not in scope

Automatic resume from the persisted vectors. The vectors make a manual restart
possible; wiring a resume path is a separate change and would need its own plan.
