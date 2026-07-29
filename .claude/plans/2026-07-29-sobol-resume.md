# Sobol: resume from persisted QoI vectors

## Goal

A Sobol run that hits its wall can be resubmitted and pick up where it stopped,
skipping evaluations already on disk. At N=256 each design matrix is 256 runs of
roughly a minute, so a resumed job that reuses `y_A` and `y_B` alone starts
20% ahead; one that reuses six of eight `y_AB<i>` starts 80% ahead.

This is what makes the study finishable at all: observed throughput on cosmos is
~9 effective cores against a 24 h wall, and N=256 does not reliably fit in one
window. Chaining windows is the only route that does not require a smaller N.

Done means: rerunning with `resume=true` against a populated `results/sobol/`
recomputes nothing that is already there, and produces indices identical to an
uninterrupted run.

## The correctness risk this must not create

Cached vectors are only valid for the design matrices that produced them. `A`
and `B` come from `MersenneTwister(seed)` and depend on `N` and `k`; the QoI
depends on `ncells`, `tpost`, `nframes`, `tmature`. Resuming across a change in
any of those silently mixes incompatible evaluations and yields indices that
look plausible and are wrong.

So resume is **guarded**, not best-effort: the run configuration is written to
`results/sobol/run_config.csv` on first write, and a resume that disagrees with
it aborts rather than proceeding.

## Steps

1. **Write `run_config.csv`** on first run — `N`, `k`, parameter names, `ncells`,
   `tpost`, `nframes`, `tmature`, `seed`.
   Verify: file appears after a fresh run and lists all eight fields.

2. **Add `resume` kwarg plus CLI flag**, default `false`. On resume, compare the
   manifest against the current configuration and `error()` on any mismatch.
   Verify: a resume with a changed `ncells` aborts with a message naming the
   field, and leaves existing files untouched.

3. **Load cached vectors when present and correctly sized**, else compute and
   persist. Applies to `y_A`, `y_B`, and each `y_AB<i>`.
   Verify: a resume against a complete directory runs zero evaluations and
   finishes in seconds.

4. **Equivalence check** — a run interrupted after 3 parameters and resumed must
   produce the same `sobol_indices.csv` as an uninterrupted run at the same seed.
   Verify: byte-identical output on an N=2 pair.

## Known limitation

The non-OK ("clog") counter only counts evaluations performed in the current
session; cached runs are not re-tallied. The printed count is therefore
per-session, not cumulative, and is labelled as such. Persisting it would mean a
second manifest write per matrix for a diagnostic number that no result depends
on.

## Files touched

- `julia/analysis/sobol_sensitivity.jl`
