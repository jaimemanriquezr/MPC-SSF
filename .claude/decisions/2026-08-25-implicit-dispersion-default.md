# Implicit dispersion becomes the default (2026-08-25)

## The decision

`parameters.ImplicitDispersion` defaults to `true`. Approved by Jaime 2026-08-25.

## Why (measured, not argued)

`probeImplicitDispersion(NCells=500, Kappa=1e-6, Days=3)`, cosmos job 3537570_2,
with `MaxDt = 1e-2` so the CFL decides rather than the cap:

| | explicit | implicit |
|---|---|---|
| steps | 1,448,486 | 203,600 |
| wall | 2289 s | 460 s |
| accuracy vs explicit | -- | 6.10e-05 |
| binding region | dispersion | enclosed L, 92% of steps |
| term shares | -- | w_v/dz 0.53, **w_a/dz^2 0.00**, w_b 0.01, w_s 0.46 |

7.11x fewer steps, 4.97x faster, for a 6.1e-05 difference -- and that figure is an
UPPER bound, since Check A's `small = 3e-6` is not CFL-free at N=500 so the two
arms did not run at identical dt.

Independent confirmation at 20 d, N=100: the implicit run reproduces the explicit
baseline to four decimals (bedPeak 0.2771, supPeak 0.3246, supFrac 79.4%,
top2/bed 10.1%). The scheme does not change the answer, only the cost.

## Alternatives rejected

- **Leave it off.** Rejected on the numbers above. At N=1000 the gap widens
  further: with dispersion explicit dt scales as dz^2 (cost ~ N^3), with it
  implicit the bound is 53% advective and 47% dz-independent so dt scales as dz
  (cost ~ N^2). Estimated 23 min vs 4.9 h for a 3 d run.
- **Believe the earlier 0.02 d result** (1.00x step reduction, 24% slower,
  accuracyRelDiff 1.42). Rejected as a startup artefact: `L_e` in the enclosed
  liquid region blows up while those liquids are near zero, so `w_s` spuriously
  dominates and dispersion cannot be the binding term. `slurm/impdisp.sbatch`'s
  header had already made this argument; it was written and never run.

## Known divergence, accepted

`SSF.jl` has NO implicit dispersion. `SSF.jl/test/golden/export_reference.m` calls
`simulate()` without the flag and therefore inherits the new default, so the
MATLAB golden references move and Julia cannot follow. CLAUDE.md names
`SSF.jl/test/golden/` the source of truth for cross-implementation agreement, so
this is a real cost, taken deliberately. Either the Julia side gains the option or
the golden export pins `ImplicitDispersion=false` explicitly.

## Caveat carried forward

This does NOT license raising `AdaptiveMaxDt`. A separate experiment the same day
raised the cap 3e-6 -> 1e-5 on the strength of a 3 d accuracy test (<=1.13%) and
the 20 d arms then aborted with Flag = "BIOFILM" at t = 15.5 d. Accuracy and
stability were tested in different regimes. The cap stands at its previous value;
see the comment at `simulate.m` `AdaptiveMaxDt`.
