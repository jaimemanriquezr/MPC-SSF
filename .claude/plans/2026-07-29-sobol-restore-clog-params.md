# Sobol: restore velocity and beta_porosity to the design

## Goal

`sobol_sensitivity.jl` estimates indices for eight parameters and excludes three —
`velocity`, `beta_porosity`, `zeta_0` — on the grounds that their ranges drive the
filter into a clogged state and leave the QoI undefined.

That was true under the old cohesion law. It is not true now: the corrected
potential removed clogging entirely, and `log_oat_{pulse,startup,flowstep}/measures.csv`
report `flag = OK` for all 27 parameters, the three excluded ones included.

Meanwhile the corrected Log-OAT ranks `beta_porosity` **first** in `pulse`
($I_"rms"$ 0.566) and `velocity` third (0.515). So the variance-based study
currently omits ranks 1 and 3 of its own top four and cannot address the central
claim of the co-author brief.

Done means: `k = 10`, with `velocity` and `beta_porosity` sampled over the same
ranges Morris uses, and a completed run whose indices are directly comparable to
the Log-OAT ranking.

`zeta_0` stays out. It ranks 0.045 in `pulse`, an order below the others, and
adding it costs a design column for a near-zero result.

## Ranges — deliberately matched to Morris

| parameter | lo | hi | scale | source |
|---|---|---|---|---|
| `velocity` | 2.16 | 21.6 | log | Schijven2013 T1 filtration-rate spread, nom 7.2 m/d |
| `beta_porosity` | 0.95 | 0.995 | lin | Diehl2025/Lund biofilm porosity, nom 0.99 |

Matching `morris_screening.jl` exactly is the point: the two studies then sample
the same distribution and their rankings can be compared without a
reparameterisation argument.

**`beta_porosity` is sampled as β itself, not as the gap.** The Log-OAT helper
`set_beta_gap()` takes `1 − β`, and log-uniform sampling of the gap is a
different distribution over β than uniform sampling of β. Since Sobol indices
depend on the input distribution, reusing the gap helper would silently make the
Sobol and Morris results incomparable. A direct setter is used instead.

## Steps

1. **Add a direct β setter** and the two `SP` entries; correct the stale
   exclusion rationale in the header comment.
   *Verify*: `k == 10`, and a fresh N=2 run completes with both new parameters
   appearing in the output table.

2. **Confirm neither new parameter clogs at its range extremes** — the original
   exclusion existed for a reason and the fix is what removed it.
   *Verify*: evaluations at the lo and hi corners of both parameters return
   `flag = OK`.

3. **Cancel the two in-flight Sobol jobs** (`3427791`, `3428726`). Both compute
   the eight-parameter design and are now obsolete.
   *Verify*: `squeue` shows neither.

4. **Resubmit at k=10**, 16 threads, 48 h wall, resume enabled.
   *Verify*: the job writes `y_A.csv` and flushed progress within the first hour,
   which the previous run never did.

## Cost

`N(k+2)` grows slowly in `k`: 3072 runs at k=10 against 2560 at k=8, so restoring
both parameters costs about 20% more, not 25% per parameter. This is the main
reason to do it now rather than as a second study.

## Files touched

- `julia/analysis/sobol_sensitivity.jl`
- `slurm/sobol-resume.sbatch` (cosmos, and see the separate task to track it in git)
