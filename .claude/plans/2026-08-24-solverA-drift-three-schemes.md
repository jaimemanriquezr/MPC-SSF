# Solver A drift: arm the probe for all three cohesion schemes

## Goal

Decide whether `matched` still earns its place alongside `shin` and `bailo`.

The claim to test: since the 2026-08-24 `v_b` fix, `bailo` already has the
operator consistency that `matched` was built to provide. Bailo's u-row flux is
`M(u_i,u_j) = zeta_0 (u_up)^+ (1-u_down)^+` times `v`; Solver B now applies
`u_up * zeta_0 (1-u_down) * grad(mu)` (simulate.m:764-767) — the same operator.
If bailo's free-running drift is <= matched's, matched has no unique virtue and
the scheme set collapses to two, which removes the motivation for a scheme class
hierarchy.

Done = a three-scheme drift comparison exists as a committed figure + per-scheme
`.mat`, and the matched/bailo question is answered from measurement.

## Steps

1. **Fix the free-running diagnostic for `matched`.**
   `simulate.m:725` gated on `if ~useBailo`, so `matched` fell through to
   `lhsPar`, built from the SHIN stencil at line 718. A matched arm would have
   measured shin.
   *Verify:* `checkcode` clean; the matched arm's free state must differ from
   shin's in the results (if `matched` and `shin` free-running curves coincide
   exactly, the branch is not being taken).

2. **Split the probe into runner + plotter.**
   New `analysis/probes/probeSolverADrift.m` runs ONE scheme and writes
   `analysis/probes/data/solverA_drift_<scheme>_n<N>.mat`.
   `analysis/plotSolverADrift.m` becomes the plot half, loading those caches and
   generalised from 2 hardcoded schemes to N.
   *Verify:* `checkcode` clean on both; plotter reproduces the existing figure
   layout.

3. **Fix the peaks-overlay bug in the plotter.**
   `D = D{k}` inside `for k = 1:2` clobbered the cell array on the first
   iteration, so the second threw. The peaks overlay has never rendered, and the
   `cellfun(@max, D)` in the last tile would also have failed. Renamed to `Dk`.
   *Verify:* the figure now carries peak markers on both/all curves.

4. **Run the three schemes as a Slurm array on cosmos.**
   `slurm/solverA_drift.sbatch`, `--array=0-2`, one scheme per task, so wall
   clock is the slowest scheme not their sum. N=100, kappa=1e-6, 20 d — matching
   the existing `drift20d.sbatch` evidence base.
   *Verify:* three `.mat` files return; all three `Flag == "OK"`.

5. **Plot locally and read off the answer.**
   *Verify:* the one-step and accumulated drift columns rank the three schemes.

## Files touched

- `src/@State/simulate.m` — matched branch in the free-running diagnostic
  (inside `if parameters.RecordCflBudget`, so diagnostic-only)
- `analysis/probes/probeSolverADrift.m` — new
- `analysis/plotSolverADrift.m` — rewritten as the plot half
- `slurm/solverA_drift.sbatch` — new

## Open

Days=20 chosen to match `drift20d.sbatch`. The 3 d numbers already show the
drift peaks creeping upward rather than saturating, so a short run would not
separate the schemes on the quantity that matters.
