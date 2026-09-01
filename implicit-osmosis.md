# Implicit integration of the osmosis term

Notes on the `implicit_osmosis` option in `julia/src/simulate.jl`. Default is `false`
(explicit forward Euler, bit-for-bit with the MATLAB path); setting it `true` integrates the
osmotic relaxation with backward Euler.

## The term

Osmosis is a single stiff linear relaxation source on the enclosed-water volume fraction `φW`,
driving the enclosed fraction `φe` toward osmotic equilibrium `φe = β·φB`:

```
rhsEnclWater = (β·φB − φe) / τ          # simulate.jl
```

with `τ = model.osmosis_rate` the relaxation time and `β = model.biofilm_porosity`.

## Affine structure in φW

Both `φB` and `φe` contain `φW` additively:

```
φB = φW + (matrixP + enclosedP)/ρP + enclosedL/ρL
φe = φW + enclosedP/ρP + enclosedL/ρL
```

Collecting the `φW`-dependence of the rate, the `φW` coefficients are `β` (from `β·φB`) and `−1`
(from `−φe`), so the source is affine in `φW`:

```
rhsEnclWater = A_osm − k_osm · φW,     k_osm = (1 − β)/τ,
```

where `A_osm` gathers the remaining (matrix/particulate/liquid) contributions, treated frozen
over the step. `k_osm > 0` is the stiff decay rate; since `β` is close to 1, `k_osm` is smaller
than the raw `1/τ`, but still the stiffest scale in the problem.

## Backward-Euler update

For the linear part `dφW/dt = A_osm − k_osm·φW` (transport handled separately and explicitly),
backward Euler gives

```
φW^{n+1} = (φW^n + Δt·A_osm) / (1 + Δt·k_osm)
        = φW^n + Δt · (A_osm − k_osm·φW^n) / (1 + Δt·k_osm).
```

The second form shows the update is exactly the explicit source **damped by `1/(1+Δt·k_osm)`**,
because `rhsEnclWater` evaluated at the current state already equals `A_osm − k_osm·φW^n`. That
is precisely how it is implemented — the whole RHS is divided by the damping factor:

```
kosm          = (1 - beta) / tau
rhsEnclWaterA = implicit_osmosis ? rhsEnclWater ./ (1 .+ dt .* kosm) : rhsEnclWater   # Solver A RHS
rhsEnclWaterB = implicit_osmosis ? rhsEnclWater ./ (1 .+ dt .* kosm) : rhsEnclWater   # Solver B φW update
```

Then the `φW` update proceeds with the damped source added explicitly:

```
phiW = phiW + (dt/dz)·(fWatIn − fWatOut)/porosity_centers + dt·rhsEnclWaterB
```

The transport (advective) flux stays explicit; only the stiff reaction source is made implicit
(an IMEX split of the `φW` equation). Applying the same damping to `rhsEnclWaterA` keeps the
Solver A biofilm-volume RHS consistent with the damped relaxation.

## Properties

- **Unconditionally stable** in the relaxation: `1/(1+Δt·k_osm) ∈ (0,1)` for any `Δt > 0`.
- **Exact as `Δt → 0`** (damping factor → 1, recovering the explicit source).
- **Equilibrium-preserving**: at `φe = β·φB` the source is zero, and damping zero is zero, so
  steady states are unchanged; `τ` and `β` retain their meaning.
- **CFL relief**: the explicit branch carries `1/τ` in the enclosed-liquid CFL weight,
  `w_b[3] = max(transport_liquid_rate/β, 1/τ)`; the implicit branch drops the `1/τ`,
  `w_b[3] = transport_liquid_rate/β`, so the step is no longer osmosis-bound. With `τ` in the
  range `1e-7…1e-3` across the presets, this is typically the binding constraint removed.

## Caveat

Solver A nominally runs on the previous step `Δt_A = Δt^n` (split-step convention), but the
`rhsEnclWaterA` damping uses the current `dt`. Since Solver A only uses this to form the
biofilm-volume RHS, the effect is minor; align the damping factor with `Δt_A` if exact
split-step semantics are required.
