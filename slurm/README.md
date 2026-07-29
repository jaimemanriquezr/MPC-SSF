# Slurm job definitions

These are the cosmos (LUNARC) batch scripts. **This directory is the source of
truth**; the copies under `~/MPC-SSF/slurm/` and `~/MPC-SSF-sobol/slurm/` on
cosmos are rsync targets, not a second home.

Until 2026-07-29 these files existed *only* on cosmos and were tracked nowhere.
Every job definition — including the ones that produced committed results — was
one `rm` or one reimage away from being unrecoverable, and two copies of
`sobol-resume.sbatch` had already diverged without anything noticing.

## Deploying

```sh
rsync -av -e "ssh -o ControlPath=$HOME/.ssh/cm-cosmos" \
      slurm/ cosmos:MPC-SSF/slurm/ --exclude logs
```

`~/MPC-SSF` on cosmos is **not** a git repo — it is an rsync target, so edit here
and push, never the other way round. See the `cosmos-hpc-access` note for how the
multiplexed connection works; Claude cannot authenticate on its own.

## What each script runs

| script | job | notes |
|---|---|---|
| `sobol.sbatch` | variance-based indices | 48 threads, 24 h. **Superseded** — see below |
| `sobol-resume.sbatch` | same, patched | 16 threads, 48 h, resumable. Paths point at `~/MPC-SSF-sobol` |
| `sobolcost.sbatch` | per-evaluation cost probe | the measurement behind the thread-count choice |
| `logoat.sbatch` | logarithmic OAT, 3 scenarios | array job |
| `lightsweep.sbatch` | irradiance amplitude sweep | array 1-5, 500 cells, `dt = 3e-7` |
| `lightprobe.sbatch`, `lightdiag.sbatch` | light submodel probes | |
| `pkgtest.sbatch` | `Pkg.test` on the cluster | |

## Two lessons these scripts encode

**Fewer threads, not more.** `sobol.sbatch` asked for 48 and sustained a
`CPULoad` near 9 — the surplus threads contend in garbage collection. The cost
probe put the whole task at 36.9 core-hours; the 48-thread job burned 147 and
finished nothing. `sobol-resume.sbatch` asks for 16.

**`--exclusive` queues badly.** With 147 of 186 `lu48` nodes allocated, an
exclusive 48-core request was scheduled 21 h out while the same job at `-c 16`
non-exclusive started the same day. `lu48` also allows 7 days, so a 24 h wall is
a choice rather than a limit — prefer a longer wall over a job that dies at one.
