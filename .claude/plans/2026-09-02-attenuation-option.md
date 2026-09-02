# Thread an `Attenuation` option through modelLund -> pathogenModel -> probeChain

## Goal

`attenuationParticle` (nu_P) is hard-coded in `src/presets/modelLund.m`. That makes it
impossible to run a control arm at the old value under today's code, so the effect of
today's nu_P change (0.094 -> 52, commit b1fcd03) cannot be separated from the effect of
`mu_HET` (2.0 -> 4.8) or from the Sept-1 kinetics/solver audit (`ee0e4d3`, 523 lines of
simulate.m, committed AFTER the chain_fld2x_* legs were grown on Aug 26 but uncommitted
at the time, so git cannot date it).

Done means: `probeChain(..., Attenuation=0.094)` runs today's code with the old
coefficient, every existing caller is byte-identical in behaviour, and a 2x2 over
{mu_HET 2.0, 4.8} x {nu_P 0.094, 52} can be submitted.

## Steps

1. `src/presets/modelLund.m` — add `options.Attenuation (1,1) {mustBeNumeric} = 52;` to
   the arguments block and use it in place of the literal for `attenuationParticle`.
   Default MUST be 52, the current preset value, so existing callers are unchanged.
   Verify: `grep -n "attenuationParticle = options.Attenuation"`.
2. `src/presets/pathogenModel.m` — add `options.Attenuation (1,1) {mustBeNumeric} = NaN;`
   (NaN = preset, the convention probeChain already uses for KDOM/KHPO4/MuHET) and
   forward it to the `modelLund(...)` call only when it is not NaN.
   Verify: `pathogenModel()` and `pathogenModel(Attenuation=0.094)` both construct.
3. `analysis/probes/probeChain.m` — add `opts.Attenuation (1,1) double = NaN;` and pass
   it to `pathogenModel`. Verify: the option appears in the leg-1 record guard.
4. `slurm/tenore_2x2.sbatch` — array 0..3, lit, 30 d, one arm per (mu_HET, nu_P) pair.
   Verify: four distinct tags, no pre-existing legs under any of them.

## Files touched

- `src/presets/modelLund.m`
- `src/presets/pathogenModel.m`
- `analysis/probes/probeChain.m`
- `slurm/tenore_2x2.sbatch`

## Verification that this changes nothing by itself

`pathogenModel()` with no arguments must produce a model identical to the one before this
plan. The golden suites in SSF.jl are the cross-implementation gate but do not cover the
MATLAB presets; the practical check is that the defaults are the current literals (52 in
modelLund, NaN passthrough elsewhere) and that no call site passes the new option.

## Supersedes

Cosmos 3564903 (`slurm/tenore_chain.sbatch`), submitted before this confound was found.
Its lit arm duplicates 2x2 arm 3, and its tags would leave a near-duplicate chain on disk
under a different name. Cancel it.
