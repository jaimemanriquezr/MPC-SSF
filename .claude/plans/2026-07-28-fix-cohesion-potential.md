# Plan: fix the Cahn-Hilliard potential to the published form

## Goal

`CahnHilliardModel`'s default `potential_gradient` does not implement the model published in the
manuscript, and never has. Fix it, prove the fix against the legacy reference, and establish a
golden case that can actually detect the difference (none currently can).

**Done** means: the Julia default equals the legacy/published form exactly; `zeta_1` becomes a live
parameter; a new golden discriminates the two forms; and the blast radius on existing goldens and on
already-published sensitivity results is measured and reported.

## The defect

Manuscript `model.tex:144` defines `Ψ = ½φ³(φ/2 − ζ₁)`, so `Ψ′ = φ²(φ − 1.5·ζ₁)`.
All three legacy model files implement exactly that, baking ζ₁ into a stored handle:

| file | particles | kappa | zeta_0 | zeta_1 | `gradient_potential` |
|---|---|---|---|---|---|
| `SDparameters.mat` | 3 | 1e-7 | 1000 | 0.005 | `@(phi) phi.^2.*(phi-3*(1/200)/2)` |
| `SDparameters_pathogen.mat` | 4 | 1e-7 | 1000 | 0.005 | same |
| `thesis_model.mat` | 4 | 1e-6 | 1 | 0.01 | `@(phi) phi.^2.*(phi-3/2*(1-0.99))` |

MPC-SSF's Julia uses the `CahnHilliardModel` default
`0.25·u²(1−u)²` (`julia/src/cohesion/CahnHilliardModel.jl:28`) and **no preset or analysis script
anywhere overrides it** (verified by grep across `src/presets/` and `analysis/`). The MATLAB twin
carries the same wrong default (`src/cohesion/@CahnHilliardModel/CahnHilliardModel.m:7`).

Measured discrepancy — not a small one, and it changes sign:

| φ | published/legacy | current default | ratio |
|---|---|---|---|
| 0.01 | 2.50e-7 | 2.45e-5 | **98×  too large** |
| 0.05 | 1.06e-4 | 5.64e-4 | 5.3× too large |
| 0.10 | 9.25e-4 | 2.03e-3 | 2.2× too large |
| 0.30 | 2.63e-2 | 1.10e-2 | 0.42× |
| 0.63 | 2.47e-1 | 1.36e-2 | **18× too small** |

φ=0.63 is the documented self-limiting peak (`worklog-2026-07-20.md:101`), so cohesive resistance to
compaction is ~18× too weak exactly where it matters most. Crossover is near φ≈0.2.

## Why no test caught it

`reference/` (SimpleModel) and `reference_adaptive/` (modelLund) both start from a **clean initial
condition**, so φ_b = 0, and *both* forms vanish at φ=0 — they cannot discriminate. The pathogen
goldens do seed φ=0.05 but run only `SimulationTime=1e-5` d (~0.86 s simulated, 10 steps), far too
short for a cohesion difference to propagate above the 1e-7 tolerance. So the blind spot is
structural, not accidental.

## Steps

1. Change the default to `u -> u.^2 .* (u .- 1.5*zeta_1)` in
   `julia/src/cohesion/CahnHilliardModel.jl`. `Base.@kwdef` permits a later default to reference an
   earlier field — verified. Update the docstring and the header comment, which currently document
   the old form.
   *Verify*: `CahnHilliardModel(kappa=1e-7, zeta_0=1e3, zeta_1=0.005).potential_gradient(φ)` equals
   the legacy handle at φ ∈ {0.01, 0.05, 0.1, 0.3, 0.63} to machine precision. (Already confirmed
   numerically identical in a standalone check.)

2. Mirror it in `src/cohesion/@CahnHilliardModel/CahnHilliardModel.m` (both the property default and
   the `arguments` default). MATLAB cannot reference another argument in a default, so build the
   handle in the constructor body after `Zeta1` is known.
   *Verify*: same five-point comparison in MATLAB.

3. **Do not** silently change `zeta_1`'s own default. It is currently `0.0`, which under the new form
   gives `Ψ′ = φ³` — a legitimate degenerate case, and the right behaviour for a caller who does not
   supply ζ₁. Presets set it explicitly (`modelPathogen.jl:133` = 1/100, `modelLund.jl:87` = 1/100).
   Note the publication models use **0.005**, not 0.01 — flag the mismatch, do not "fix" it unasked.

4. Fix `set_zeta0()` (`analysis/log_oat_sensitivity.jl:111`) to carry `mobility` and
   `potential_gradient` through the reconstruction instead of dropping them. This is currently
   harmless-by-accident; after step 1 it becomes load-bearing, because dropping
   `potential_gradient` is exactly what must happen for `set_zeta1()` to work and exactly what must
   *not* happen for `set_zeta0()`. Add `set_zeta1()` that deliberately omits it so the handle
   regenerates from the new ζ₁, and `set_kappa()` (a ready-made `b_kappa()` sits as dead code at
   `analysis/morris_screening.jl:102-103`).
   *Verify*: after `set_zeta1(m, 0.02)`, `m.cohesion_submodel.potential_gradient(0.05)` reflects
   ζ₁=0.02, not the old value. This is the trap that would silently produce a zero sensitivity.

5. Add a discriminating golden: seed a mature biofilm and run long enough for cohesion to act. Export
   from the legacy code with `export_pathogen_reference.m` (the existing cross-repo driver) at a
   longer `SimulationTime`, and confirm the new reference **fails** against the old default and
   **passes** against the fixed one. A test that cannot fail against the bug is not a test.

6. Measure the blast radius and report before anything is published:
   - which existing goldens move, and by how much;
   - re-run one representative log-OAT scenario and diff the ranking against the committed
     `measures.csv`. Every Julia sensitivity result to date — the log-OAT rankings, the Sobol
     indices, the `zeta_0` clog-driver classification — was computed with the wrong cohesion law.
     Attachment dominance plausibly survives (attachment is independent of cohesion) but the biofilm
     block almost certainly shifts.

## Files touched

- `julia/src/cohesion/CahnHilliardModel.jl` — the default + docs
- `src/cohesion/@CahnHilliardModel/CahnHilliardModel.m` — MATLAB twin
- `julia/analysis/log_oat_sensitivity.jl` — `set_zeta0` fix, new `set_zeta1`/`set_kappa`
- `julia/test/golden/` — new discriminating reference
- `.claude/decisions/2026-07-28-published-cohesion-potential.md` — the why

## Risk

This changes simulated physics. Any manuscript figure produced by the **Julia** implementation rather
than the legacy code may move; which implementation produced which figure must be checked before any
claim is made. Runs expected to exceed ~1 minute go to cosmos per standing instruction.
