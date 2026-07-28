# The Cahn-Hilliard potential default is the published Ψ′, not the double well

## The decision

`CahnHilliardModel`'s default `potential_gradient` is now

```julia
u -> u .^ 2 .* (u .- 1.5 * zeta_1)
```

replacing `u -> 0.25 .* (u .^ 2 .* (1 .- u) .^ 2)`. `zeta_1` becomes a live parameter for the first
time. The MATLAB twin (`src/cohesion/@CahnHilliardModel/CahnHilliardModel.m:7` and
`src/@Model/Model.m:36`) needs the same change — it is **still outstanding**, see "Follow-up".

## Why

The manuscript (`manuscripts/AWR-SSF/model.tex:144`) defines

```
Ψ(φ) = ½ φ³ (φ/2 − ζ₁)     ⇒     Ψ′(φ) = φ² (φ − 3ζ₁/2)
```

and `results.tex:80` tabulates ζ₁ = 1e-2 as a model parameter with a citation. **All three legacy
model files implement exactly that**, baking ζ₁ into a stored function handle rather than reading the
`zeta_1` field (verified by `func2str` under MATLAB R2025b):

| file | particles | kappa | zeta_0 | zeta_1 | `gradient_potential` |
|---|---|---|---|---|---|
| `SDparameters.mat` | 3 | 1e-7 | 1000 | 0.005 | `@(phi) phi.^2.*(phi-3*(1/200)/2)` |
| `SDparameters_pathogen.mat` | 4 | 1e-7 | 1000 | 0.005 | same |
| `thesis_model.mat` | 4 | 1e-6 | 1 | 0.01 | `@(phi) phi.^2.*(phi-3/2*(1-0.99))` |

The Julia port never overrode `potential_gradient` — not in any preset, not in any analysis script —
so it silently used a **different functional form** with no ζ₁ and different equilibria (0 and 1,
rather than 0 and 3ζ₁/2). The error is large and changes sign across the physical range:

| φ | published | old default | ratio |
|---|---|---|---|
| 0.01 | 2.50e-7 | 2.45e-5 | **98× too large** |
| 0.05 | 1.06e-4 | 5.64e-4 | 5.3× too large |
| 0.10 | 9.25e-4 | 2.03e-3 | 2.2× too large |
| 0.30 | 2.63e-2 | 1.10e-2 | 0.42× |
| 0.63 | 2.47e-1 | 1.36e-2 | **18× too small** |

Crossover is near φ ≈ 0.2. φ = 0.63 is the documented self-limiting peak
(`worklog-2026-07-20.md:101`), so cohesive resistance to compaction was ~18× too weak exactly where
the physics matters most. This is not a tolerance-level discrepancy; it is a different cohesion law.

After the fix, the default reproduces both legacy handles to machine precision — worst absolute
difference **0.0** against `SDparameters*.mat` (ζ₁=0.005) and **2.2e-19** against `thesis_model.mat`
(ζ₁=0.01, float rounding on the `3/2*(1-0.99)` literal). With ζ₁=0 it degenerates to `u³`, which is
the correct behaviour for a caller who supplies no ζ₁, and is why `zeta_1`'s own default was left at
`0.0` rather than being quietly set to a preset value.

## Why no existing test caught it

The blind spot is structural, not accidental:

- `reference/` (SimpleModel) and `reference_adaptive/` (modelLund) both start from a **clean initial
  condition**, so φ_b = 0 — and *both* forms vanish at φ = 0. They cannot discriminate at any
  tolerance.
- `reference_pathogen*/` does seed a mature biofilm (φ = 0.05) but runs `SimulationTime = 1e-5` d
  (~0.86 s simulated, 10 steps), far too short for a 5× difference in the chemical potential to
  propagate into concentrations above the `atol=1e-7` gate.

So four golden suites passed against the legacy code while implementing different physics. **A test
that cannot fail against the bug is not a test for the bug** — hence the new discriminating golden
required by the plan (seeded biofilm, horizon long enough to matter, must fail against the old
default).

## The closure trap

Because the default closes over `zeta_1`, reconstructing the struct and passing
`potential_gradient` through carries the **old** ζ₁ with it. A `set_zeta1` written the obvious way
would change the field and leave the physics untouched — reporting a clean zero sensitivity that
looks like a finding. Guarded explicitly:

- `set_kappa` / `set_zeta0` (`analysis/log_oat_sensitivity.jl`) carry `mobility` and
  `potential_gradient` through via the `_ch` helper.
- `set_zeta1` deliberately **omits** `potential_gradient` so it regenerates from the new value.
- Asserted in a standalone check: `set_zeta1(m, 0.02)` gives `dpsi(0.05) = 5.0e-5` (matching
  `0.05²(0.05 − 1.5·0.02)`), while `set_zeta0` and `set_kappa` leave `dpsi(0.05)` at the baseline
  `8.75e-5`. Note the previous `set_zeta0` dropped both handles, which was harmless only because the
  old default ignored every field.

## Alternatives rejected

**Override `potential_gradient` in each preset instead of changing the default.** Rejected: it leaves
a wrong default as a trap for the next caller, and the presets are ports of `.mat` files whose whole
point is to carry the published parameters. The default is where the published model belongs.

**Read `zeta_1` inside the solver rather than closing over it.** Cleaner in principle — it removes the
closure trap entirely — but it changes the `potential_gradient` call signature and breaks parity with
the MATLAB `function_handle` field, which is the interface the legacy `.mat` files rely on. Rejected
for now; worth revisiting if the closure trap bites again.

**Match `thesis_model.mat`'s ζ₁ = 0.01 as the default.** Rejected: the publication models use
**0.005**. Picking either as a global default would misrepresent the other. Presets set it explicitly;
the mismatch between the test model (0.01) and the publication models (0.005) is flagged for Jaime
rather than silently resolved.

## Consequences to chase down

Every Julia sensitivity result computed so far — the log-OAT rankings across three scenarios, the
Sobol indices, and the `zeta_0` clog-driver classification — used the wrong cohesion law. Attachment
dominance (`attach_sand` 0.531, `sand_pathogen` 0.480) plausibly survives, since attachment is
independent of cohesion, but the biofilm block should be expected to move. The mesh-convergence
numbers may also shift. None of this can be asserted without re-running; the plan's step 6 covers it.

Separately: if any manuscript figure was produced by the **Julia** implementation rather than the
legacy code, it may move. Which implementation produced which figure is not yet established and must
be before any claim is made.

## Follow-up (not yet done)

The MATLAB twin still carries the old default in two places
(`src/cohesion/@CahnHilliardModel/CahnHilliardModel.m:7`, `src/@Model/Model.m:36`), and it is live:
`src/@State/simulate.m:111` reads `model.CohesionSubModel.PotentialGradient` and applies it at
`:339`. It was deliberately not edited in this pass because a worktree-isolated agent held
`matlab-claude` for the pathogen port at the time. MATLAB cannot reference a sibling argument in an
`arguments` default, so the handle must be built in the constructor body after `Zeta1` is known —
declare `input.PotentialGradient function_handle` with **no** default so `isfield` can detect
omission.
