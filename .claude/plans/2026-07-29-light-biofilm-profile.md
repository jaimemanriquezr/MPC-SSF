# Plan: simplified HET/PHO model and the light–biofilm-profile study

## Goal

Jaime's hypothesis: at higher light incidence, more biofilm forms in the upper filter (z < 0) and a
thicker cake-like phototroph structure develops.

**Done** means: a simplified two-particle model runs ≥7 simulated days at `ncells ≈ 500` with biofilm
forming, swept over irradiance amplitude, reporting (a) the location of the peak biofilm volume
fraction and (b) biofilm mass above z = 0.

## The finding that reshapes the experiment

The light factor is (`src/simulate.jl:358-360`)

```
I_eff = I(t) · exp(−η) / light_optimal      lightEffective = I_eff · exp(1 − I_eff)
lightFactor = max(minimum_light_factor, lightEffective)
```

a Steele curve peaking at `I_eff = 1`. For `modelLund`, `light_optimal = 0.01814` and peak
irradiance is `0.8`. So at the top of the supernatant, with η ≈ 0:

| quantity | value |
|---|---|
| `I_eff` at supernatant top | **44.1** |
| `lightEffective` there | 8.5e-18 |
| after the floor | 0.01 = `minimum_light_factor` |

**The entire supernatant sits at the 1% dark-respiration floor — not because it is dark, but because
it is ~44× above the photoinhibition optimum.** Phototroph growth is suppressed by *too much* light.

Where the optimum actually falls, computed on the real η profile:

| ncells | dz | peak light factor | location | cells above the floor |
|---|---|---|---|---|
| 30 | 32.8 mm | 0.314 | z = 0 | **1** |
| 500 | 2.0 mm | 0.869 | **z = +2.0 mm** | **4** |

So the phototroph growth optimum is a thin band a **few millimetres inside the bed**, where
attenuation has brought `I_eff` down to ≈1.

### Consequences

1. **Jaime's cake prediction is right, but the mechanism is inverted.** A concentrated phototroph
   band does form near the sand surface — because that is where photoinhibition *relents*, not where
   light is strongest.
2. **Raising irradiance amplitude will not increase biofilm in z < 0.** It pushes the optimum
   *deeper* (more attenuation needed to reach `I_eff` = 1) while the supernatant stays pinned at the
   floor. The hypothesis as stated predicts the wrong direction.
3. **Lowering irradiance is the intervention that populates the supernatant.** `I_eff` = 1 at η ≈ 0
   needs `I ≈ 0.0181`, roughly 44× below default.
4. **`ncells ≈ 500` is vindicated and is the minimum.** At 30 cells only one cell clears the floor
   and the peak factor is 0.314 versus a true 0.869 — the band is badly under-resolved.

So the sweep must **bracket the optimum logarithmically in both directions**, not increase
monotonically:

`I_peak ∈ {0.018, 0.05, 0.2, 0.8 (default), 3.2}` — spanning optimum-at-the-surface through
deep-inhibition.

## Geometry decision: δ = 20 mm, fixed

Jaime's correction: bare-sand attachment scales with `(1 − porosity)`, and
`computeporosity = clamp((ε₀−1)/δ·z + ε₀, ε₀, 1)` ramps porosity from 1 at z = −δ to ε₀ at z = 0.
So solid fraction — and therefore attachment — is present throughout **[−δ, 0)**, above the nominal
surface. My earlier claim that the bare-sand term vanishes in the supernatant holds only for z < −δ.

δ also sets surface optical thickness linearly: `η_sand(z→0⁻) = A(1−ε₀)·δ/2`, i.e. the roughness
layer is optically equivalent to δ/2 of packed bed. Raising δ therefore widens the attachment zone
*and* deepens attenuation, dragging the near-optimal light band up out of the bed:

| δ | light optimum | porosity there | roughness cells @500 |
|---|---|---|---|
| 5 mm (default) | +2.0 mm (in bed) | 0.40 | 2 |
| **20 mm (chosen)** | **−8.0 mm** | **0.64** | **10** |
| 100 mm | −71.9 mm | 0.83 | 50 |

**δ = 20 mm held fixed**, irradiance swept. At that geometry the near-optimal band sits 8 mm above
the nominal surface in water that is 36% solid, so phototrophs can photosynthesise near-optimally
*and* attach — the configuration the cake hypothesis requires. 20 mm is defensible for a filter
surface carrying accumulated deposit; 50–100 mm would need separate justification.

Caveat to respect: δ enters the QoI's own weighting, since `Mb = Σ(porosity·φ_b·dz)`. Holding it
fixed across the sweep keeps the comparison clean.

## The simplified model

Agreed scope: particles HET, PHO; liquids O2, DOM, IC. Built as a builder function in
`julia/analysis/` (following the `pathogen_repro.jl` precedent) rather than a new `src/presets/`
entry, to avoid a package API change and its test burden for an exploratory model.

Derived from `modelLund()` by dropping NH4, HPO4, POM and rerouting death straight to DOM:

| Reaction | Order | Monod | Stoichiometry |
|---|---|---|---|
| HET growth | HET | O2 3e-3, DOM 2e-4 | HET +1, O2 −1.2317, IC +0.3848, DOM −1.5873 |
| PHO growth (light) | PHO | IC 2e-5 | PHO +1, O2 +0.9301, IC −0.3600 |
| HET death | HET | — | HET −1, DOM +0.9123, O2 +0.0234 |
| PHO death | PHO | — | PHO −1, DOM +0.6316, O2 +0.2005 |

Retains `minimum_light_factor = 0.01`, `optimal_light_factor = 1.814e-2` and the Lund rates so the
light behaviour stays comparable. Dropping the NH4/HPO4 Monod terms makes it non-conservative in N
and P — acceptable since neither is tracked, and it removes the HPO4 starvation artifact flagged in
`SESSION_LOG_2026-07.md:165` (influent HPO4 = 0 against K = 1.4e-8). Hydrolysis goes with POM, so DOM
is recycled directly from death instead of via particulate matter.

Note this is **not** a speed optimisation. Measured legacy cost is exactly linear in cell count
(1002 / 2019 / 4042 s per simulated day at 20 / 50 / 100 cells), which proves the step size is pinned
by `adaptive_max_dt`, not by a CFL or half-saturation bound. The simplification buys interpretability,
not wall-clock.

## Steps

1. Write `julia/analysis/light_profile_model.jl` with the builder above plus an influent vector
   covering O2, IC, DOM. *Verify*: model builds, 2 particles / 3 liquids / 4 reactions; one short run
   returns `flag = OK`.
2. Measure Julia cost at `ncells = 500`. The `run_proxy` dt cap of 3e-6 was tuned for ~30 cells, so
   the CFL may genuinely bind here and the legacy linear cost model will not transfer. *Verify*: a
   0.1-day timing run gives seconds per simulated day, and the realised dt is recorded — this also
   settles the dt-limiter question at the mesh that matters.
3. Confirm biofilm actually forms and, critically, **whether φ_b in z < 0 ever becomes
   non-negligible**. `zero_cell_analysis.md` reports the supernatant contributing ~0 to total mass,
   and in the supernatant porosity = 1 so the bare-sand attachment term vanishes, leaving only
   autocatalytic biofilm-on-biofilm attachment. If no supernatant biofilm forms from a clean start,
   no light forcing can produce a cake, and that is a model-structure result. *Verify*: report
   max φ_b for z < 0 over the run.
4. Run the 5-point irradiance sweep, ≥7 simulated days, `ncells = 500`, as a cosmos sbatch array
   (one task per amplitude). *Verify*: all tasks `flag = OK`; QoIs written to CSV.
5. Report the two QoIs versus amplitude: argmax φ_b location, and ∫φ_b over z < 0. Expect the peak to
   move *deeper* with increasing amplitude — the opposite of the hypothesis — and supernatant mass to
   *fall*. If it does not, the reasoning above is wrong and that is the more interesting result.

## Files touched

- `julia/analysis/light_profile_model.jl` — new
- `julia/analysis/light_profile_sweep.jl` — new driver
- `julia/analysis/light_profile_report.jl` — new, step 5 aggregator
- `slurm/lightsweep.sbatch` on cosmos — new
- `julia/analysis/results/light_profile/` — outputs

## Step 5 result (measured 2026-07-29, job 3428150, all five tasks COMPLETED)

7 sim-days, 500 cells, δ = 20 mm, dt = 3e-7, all `flag = OK`, 4:47–5:40 wall each.

| amplitude | z_peak (mm) | φ_b peak | M(z<0) | M(total) | frac above |
|---|---|---|---|---|---|
| 0.018 | +2.00 | 0.02673 | 2.367e-4 | 3.405e-3 | 0.070 |
| 0.05 | +2.00 | **0.02679** | **2.380e-4** | 3.408e-3 | 0.070 |
| 0.2 | +2.00 | 0.02243 | 1.861e-4 | 2.764e-3 | 0.067 |
| 0.8 (default) | +2.00 | 0.02110 | 1.727e-4 | 2.580e-3 | 0.067 |
| 3.2 | +2.00 | 0.02079 | 1.631e-4 | 2.530e-3 | 0.064 |

**Jaime's hypothesis is contradicted.** Raising irradiance does not put more biofilm above
z = 0. Supernatant mass *falls* 31% from its optimum, and total mass falls 26%.

**There is a genuine interior optimum at A ≈ 0.05**, close to the A ≈ 0.018 predicted for
`I_eff` = 1 at the supernatant top. Both φ_b peak and M(z<0) are non-monotonic — 0.05
slightly exceeds 0.018 before the decline. This is the photoinhibition signature, measured.

### But QoI (a) — the location of the peak — is the wrong observable

`z_peak` is **+2.00 mm at every amplitude**, and that is not a coincidence or a
sub-grid shift. The profile around it:

```
z = -2.0 mm : 0.01493      z = +2.0 mm : 0.02673   <- argmax
z =  0.0 mm : 0.01508      z = +4.0 mm : 0.02652
                           z = +6.0 mm : 0.02633
```

φ_b nearly doubles across the sand surface (0.0151 → 0.0267, a 77% step over one cell)
and then decays *monotonically* downward. There is no interior maximum. The argmax is
simply the first in-bed cell, and it is pinned there by the porosity ramp — a geometry
feature — for every amplitude.

So `argmax φ_b` could never have been a light-sensitive diagnostic, whatever the light
model does. The step-5 aggregator correctly refused to call this the photoinhibition
prediction. **Use M(z<0), or a mass-weighted centroid, as the location QoI instead** — the
top-quartile centroid does move, 164.3 → 156.2 mm, though it moves *up* with amplitude
rather than down.

The prediction in "Consequences" item 2 above — that the peak moves deeper — is therefore
not tested by this sweep and needs a different observable to test at all.

## Step 2 result (measured 2026-07-29)

Cost at `ncells = 500`, `dt = 3e-7` fixed, measured locally: **~1090 s per
simulated day** warm (2345 s including compilation), i.e. **~2.1 h for the 7-day
run**. Cost per step is constant under fixed `dt`, so the extrapolation is linear
and safe.

**Cost is therefore not the binding risk** — even at 2–3× slower cosmos cores the
7-day run needs 4–6 h against a 12 h wall. The binding risk is the stability one
already flagged: `dt = 3e-7` is verified only to 0.3 d, and a run that trips the
biofilm negativity guard stops early with a partial `t_final`. The step-5
aggregator reports `flag` and `t_final` per amplitude for exactly this reason,
and refuses to draw a direction from runs that did not reach `OK`.

## Risks

- **No supernatant biofilm** (step 3). Would make the cake question moot; report as structural.
- **Cost at 500 cells is unmeasured in Julia.** If the CFL binds, cost could jump superlinearly.
  Measure before launching the sweep.
- **`minimum_light_factor` floor may dominate everything.** If phototrophs sit at the 0.01 floor
  across the whole domain except a 4-cell band, the light signal may be too weak to move the profile
  at all. That is itself a reportable finding about the light submodel's effective dynamic range.
