# The 0-cell and the mesh-dependency of the biofilm peak & clogging

Reference note on why the schmutzdecke peak biofilm fraction (`phibmax`) and the
clogging time do **not** converge under mesh refinement, while the removal QoIs
do. Grounds the caveat raised by the mesh-convergence study
(`results/mesh_convergence.csv`) in the model's treatment of the mixed cell at the
water/sand interface (Diehl et al. 2025).

## Background: the mixed 0-cell (Diehl2025)

The domain has supernatant water in `z < 0` and packed, saturated sand in
`0 ≤ z < B`. The two regions meet at `z = 0`, where porosity and phase content are
discontinuous, so a flux evaluated *at* `z = 0` is undefined. Diehl2025 avoids this
by making a **cell centre sit exactly on `z = 0`** (the "0-cell"), which straddles
the interface — half water (`z < 0`), half packed sand (`z > 0`). The grid builder
enforces this:

- `SandFilter.jl` `addgridpoints`: if a cell *face* would land on `z ≈ 0`, the cell
  count is bumped by one, so **no computational face sits at `z = 0`** — fluxes are
  only ever evaluated at faces `±dz/2`, inside single-phase regions.
- `gridzero(f)` returns the index of the cell whose centre is within `dz/2` of 0.

The mixed cell therefore represents the integrated water→sand mass exchange across
the interface; as `dz → 0` this cell shrinks onto `z = 0`.

## How the code represents water vs sand

Porosity is a linear ramp (`SandFilter.jl` `computeporosity`):

```
ε(z) = clamp( (ε0 − 1)/δ · z + ε0 ,  ε0 , 1 )
```

with `ε0 = sand_porosity = 0.40` and `δ = sand_roughness = 5×10⁻³ m`. So
`ε = 1` (pure water) for `z ≤ −δ`, ramps to `ε0` at `z = 0`, and stays `ε0` for
`z > 0`.

The **mass exchange** that builds the schmutzdecke is the sand-attachment term
(`simulate.jl:385`), weighted by the **sand volume fraction `(1 − ε)`**:

```julia
attFactorP = (1 .- porosity_centers) .* sand_factors  .+  porosity_centers .* phiBiofilm
attF       = attFactorP .* globalFlowingP .* attachment_rates      # flowing → matrix
```

`(1 − ε)` is `0` in water and `1 − ε0 = 0.60` in sand: flowing particles convert to
biofilm matrix exactly where sand appears.

## Diagnostic: the ramp is sub-grid, the peak is at the 0-cell

Mature-filter run (`pathogen_model`, `run_proxy`, `tmature = 1.5`) at three meshes
(`zerocell_diag.jl`). δ = 5 mm, ε0 = 0.40, domain z ∈ [−1, 1].

```
ncells   dz (mm)   dz/δ    0-cell(z=0) ε    peak φb   peak location
  20       48.8     9.8        0.40         0.158     IN the 0-cell (z=0)
  40       24.7     4.9        0.40         0.215     first sand cell (z=+0.025)
  80       12.4     2.5        0.40         0.282     first sand cell (z=+0.012)
```

Near-surface profile (ncells = 40), showing the porosity **step** and the peak:

```
 z = −0.049   ε = 1.00   1−ε = 0.00   φb = 0.000     (pure water)
 z = −0.025   ε = 1.00   1−ε = 0.00   φb = 0.091     (biofilm grown UP into supernatant)
 z = +0.000   ε = 0.40   1−ε = 0.60   φb = 0.143     ← 0-cell
 z = +0.025   ε = 0.40   1−ε = 0.60   φb = 0.215     ← PEAK (first sand cell)
 z = +0.049   ε = 0.40   1−ε = 0.60   φb = 0.157
```

Two facts stand out:

1. **The porosity ramp is entirely sub-grid** (`dz/δ = 2.5–9.8`). The transition
   never resolves as a ramp — it collapses to a **step**: the cell above the 0-cell
   is pure water (`ε = 1`), the 0-cell is full sand (`ε = ε0`, `1−ε = 0.60`). The
   0-cell's "half water" character is *not* represented in its (centre-evaluated)
   porosity; it is encoded only through the staggered placement of the centre on
   `z = 0`.
2. **The schmutzdecke peak sits at the 0-cell / first sand cell** — where the
   attachment weight `(1−ε)` jumps `0 → 0.60`. Biofilm also grows *upward* into the
   supernatant (`φb ≈ 0.09` at `z < 0`), carried by the Cahn–Hilliard biofilm
   velocity, consistent with Diehl2025 ("growing a couple of cm above into the
   supernatant").

## Why the peak & clog time don't converge (but Lmean & Mb do)

The peak grows ≈ `1/dz` (0.158 → 0.215 → 0.282 as dz halves twice). The mechanism:

> The **areal** biofilm mass deposited at the surface per unit time is roughly
> mesh-independent (it is a surface flux). But it is spread over a surface cell of
> thickness `dz`, so the **volume fraction** `φb ≈ (areal mass)/(dz·ρ)` scales like
> `1/dz`. The schmutzdecke is physically a near-singular thin surface accumulation:
> a coarse mesh smears it, a fine mesh sharpens it.

Consequences:

- **Integral QoIs converge.** `Lmean` and `Mb = ∫ φb dz` are ≲1–2 % by 30 cells,
  because the finite-volume *integral* of the 0-cell mass exchange is a proper
  conservative discretization — the very point of the mixed-cell construction.
- **Pointwise / threshold quantities do not.** `phibmax` (0.25→0.43 over 15→50
  cells in the convergence study) and the **clog time** — which triggers when a
  *single* surface cell crosses `φb > 0.99` (`simulate.jl:327`) — inherit the
  `1/dz` peak growth. They are resolution-limited by construction, not by a bug.

This is exactly the split the mesh-convergence note reported: the sensitivity
*rankings* (built on integral removal QoIs) are numerically robust, but the
*clogging* results are qualitative — which parameters push toward failure is
reliable, their exact timing is mesh-bound.

## Caveats & options

1. **To make clogging mesh-independent**, resolve the roughness length: `dz ≲ δ =
   5 mm` ⇒ `ncells ≳ 400` for this 2 m domain (cluster-scale) — so the porosity
   ramp and the surface accumulation are physically resolved. Alternatively, base
   the clog criterion on an integrated/averaged quantity rather than a pointwise
   cell `φb`.
2. **The peak width is also set by the Cahn–Hilliard interface parameter `κ`**,
   which couples to the mesh (diffuse-interface width). This is the UQ plan's
   warning to check "whether changing `κ` while refining the mesh changes physical
   outputs or merely the numerical interface width." It is why the clog-driving
   parameters `beta_porosity` and `zeta_0` (which enter the biofilm mechanics) were
   reported *separately* from the log-sensitivity ranking rather than mixed into it.

## Sensitivity-study implication (verified)

The right response to the sub-grid scales is **not** to restrict the whole study to
δ/κ-resolving meshes (dz ≲ 1–5 mm ⇒ ncells ≳ 400–2000, cluster-scale), but to use
only **mesh-convergent QoIs** on the feasible mesh:

- **Removal sensitivity (integral QoIs) — coarse mesh OK, verified.** A mesh-
  robustness check of the Log-OAT `I_rms` (`robustness_mesh.jl`, 30 vs 60 cells)
  shows the *sensitivity* — not just the baseline — is stable for the dominant
  parameters: sand_pathogen ×0.99, attach_sand ×1.01, influent_PAT ×1.02
  (ratios 60/30). Only the sub-0.03 tail wobbles (dispersivity ×1.29, temperature
  ×0.86), which is the numerical-noise floor the convergence study already flagged;
  the ordering is unchanged. So the removal ranking stands on the coarse mesh even
  though δ and √κ are sub-grid, because Lmean is a conservative integral and its
  sensitivity to the (surface-acting) attachment parameters is mesh-stable.
- **Clogging / peak (pointwise-threshold QoIs) — coarse mesh gives qualitative
  results only.** The clog-driving parameters are reported (which perturbations
  push to failure), but the clog *times* and `phibmax` are mesh-artifacts and must
  not be read quantitatively. A quantitative clogging study needs a ncells ≳ 400
  cluster run that resolves δ (and ideally √κ).

Rule of thumb adopted: **rank on integral QoIs; treat clogging qualitatively.**

## Pointers

- Code: `SandFilter.jl` (`computeporosity`, `addgridpoints`, `gridzero`);
  `simulate.jl:385` (sand-attachment `(1−ε)` weight), `:327` (clog trigger).
- Data: `results/mesh_convergence.csv`; diagnostic `zerocell_diag.jl`.
- Related: `SENSITIVITY_ANALYSIS.md`, Diehl et al. (2025) §2 (model), Fig. 1.
