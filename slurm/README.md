# Slurm job definitions (MATLAB)

The cosmos (LUNARC) batch scripts for the MATLAB model. **This directory is the
source of truth**; `~/MPC-SSF-manuscript/slurm/` on cosmos is an rsync target,
not a second home.

Until 2026-08-22 this directory also held the Julia campaign scripts. They moved
to `code/1d/SSF.jl/slurm/` when the Julia implementation was consolidated there
and MPC-SSF became MATLAB-only; the operational lessons about thread counts and
`--exclusive` live in that README.

## Deploying

```sh
rsync -av -e "ssh -o ControlPath=$HOME/.ssh/cm-cosmos" \
      --exclude logs . cosmos:MPC-SSF-manuscript/
```

`~/MPC-SSF-manuscript` on cosmos is **not** a git repo — it is an rsync target,
so edit here and push, never the other way round. See the `cosmos-hpc-access`
note for how the multiplexed connection works; Claude cannot authenticate on its
own. Slurm starts a job in the submit directory, so `sbatch` from the repo root:
every script here uses paths relative to it.

Probe and sweep runs go to cosmos, never to a local MATLAB batch.

## What each script runs

| script | job | notes |
|---|---|---|
| `manuscript_summer.sbatch` | E3/E1-summer/X1 at 500 cells | 90 d, 24 h wall. Builds the `mature30_summer` cache every other experiment consumes; submit stage 2 with `--dependency=afterok:<this>` |
| `manuscript_rest.sbatch` | the remaining manuscript suite | consumes that cache |
| `theta_sweep.sbatch` | `theta_growth,PHO` diagnostic sweep | array 0-3, seasons x theta |
| `mass_check.sbatch` | zero-biology mass bookkeeping + a theta asymptote | array 0-1 |

## Wall-clock note

The 500-cell 90-day baselines (jobs 3524843 / 3524845) hit the 24 h limit with
**no output written** — a 90-day run at that resolution does not fit
`--time=1-00:00:00`. `lu48` allows 7 days, so the 24 h wall in
`manuscript_summer.sbatch` is a choice rather than a limit. Either raise it or
split the run into chained 30-day stages via the mature-cache mechanism. At 100
cells the same run takes ~75 min.
