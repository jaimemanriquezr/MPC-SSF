# PLANS — backburner register

Work that is scoped and understood but deliberately parked. Distinct from `.claude/plans/`,
which holds plans for work being executed now, and from `TODO.md`, which holds small items.

Each entry states what it is, why it is parked, what unblocks it, and what was already
established so the investigation does not have to be redone.

---

## Bailo2023 implicit finite-volume Cahn–Hilliard scheme

**Status:** parked 2026-08-23. Blocked on the κ fix
(`.claude/plans/2026-08-23-recover-kappa.md`).

`documents/Bailo2023.pdf` — Bailo, Carrillo, Kalliadasis & Perez, *Unconditional
Bound-Preserving and Energy-Dissipating Finite-Volume Schemes for the Cahn–Hilliard
Equation* (arXiv:2105.05351v2). Upwinded degenerate mobility plus semi-implicit convex
splitting of the free energy. The PDF is in `documents/` and is **cited nowhere** in the
workspace.

**Why it is worth having, stated honestly.** Not for speed: the CH term is solved implicitly
and, once κ is live, will restrict dt only through `vbmax`, not through a dz⁴ stability
condition. The value is:

- **Unconditional bound preservation** `0 ≤ φ_b ≤ 1`. Nothing enforces this today —
  `src/@State/simulate.m:555-578` *aborts* with `Flag = "BIOFILM"` on negativity, and the
  golden suites avoid the Lund preset precisely because it "trips a negativity guard
  immediately" (`SSF.jl/test/golden/README.md`). Bound preservation would turn a numerical
  failure mode into a physical event, and let the goldens cover the manuscript's own ecology.
- **A discrete energy law.** The current scheme treats Ψ′ fully explicitly and has none.

**Agreed scope: standalone prototype first**, no integration.
1. Reproduce the paper's own 1-D test cases, validating the implementation independently.
2. Then substitute this model's `M(φ) = φ(1−φ)` and `Ψ'(φ) = φ²(φ − 1.5ζ₁)`, ζ₁ = 0.01.
   **Check the energy scaling before trusting the dissipation test** — this double well is
   very asymmetric, with a spinodal band only on `0 < φ < 0.01`, unlike the symmetric quartic
   the paper assumes.
3. Report: monotone energy decay? bounds preserved at large dt? Newton cost per step versus
   the current single linear solve?

**Reuse — there is an in-house reference.** `code/2d/ssf-ch-stokes-fenics` already implements
the same family: Eyre convex splitting (`main.py`, `dpsi_ey = PSI_0*(dpsi_i(u) + dpsi_e(u_0))`),
mobility upwinding around `u_m = ρ_P/(2+γ)` (`main.py:190-240`), mass lumping, and
`tests/test_potentials.py` enforcing that the split sums back to `f`. Its dolfinx port reports
"monotone energy decay, machine-precision mass conservation". Read that before writing
anything.

**Bring forward if:** the dt cost measured after the κ fix is worse than ~2×, or the
negativity abort starts tripping in production runs.

---

## Graded (non-uniform) mesh around z = 0

**Status:** parked 2026-08-23. Blocked on the κ fix — refining around an interfacial layer
is pointless while no interfacial layer is computed.

**Agreed scope: convergence studies only**, not production.

**The good news.** The physics core is already mesh-agnostic and needs no change:
conservative flux-difference form, face-weighted porosity, analytic ε (`computePorosity`)
and analytic η (`LightAttenuationEta*`), all four boundary conditions, all plotting.

**Break surface** (from the 2026-08-23 audit):

| Class | Sites | Notes |
|---|---|---|
| Silent failure | ~8 | `GridSize` telescoping mean and `GridZero` vector (both fixed by the κ work); CFL `simulate.m:491`; and four **unmarked** centre↔face averages `0.5*(x(2:end)+x(1:end-1))` at `simulate.m:150, 391, 423, 438` — these will **not** appear in a `dz` grep and become first-order on a graded mesh |
| Real work | ~10 | CH stencil assembly (shared-magnitude 4-entry stencils), `scrapeState`, `probeCover` `mat_depth`, and three gradients at `simulate.m:435, 508, 512` that need **centre spacings**, not cell widths — using `dzCell` there is a silent O(1) error |
| Mechanical | ~50 | three divergence updates `(dt/dz)→(dt./dzCell)`, `etaParticles` cumsum at `:338`, and ~45 `sum(...)*dz` quadratures across `analysis/`; `totalBiomass`'s signature changes |

**Approach.** Add `addGradedGridPoints` alongside `addGridPoints`, leaving the latter
bit-identical so all ~30 existing call sites keep reproducing published numbers. Invariants
the generator must hold: a cell **centre** exactly at z = 0 (the `+1/2` guarantees this today
and `GridZero` depends on it), the ghost cell below `Depth`, and `−Height` hit exactly.
Refine only `z ∈ [−δ, 0]` — the CH system covers cells `1..n0`, since `v_b ≡ 0` in the bed.

Also fix `etaParticles` while there: it is already first-order-inaccurate (whole cell widths
at centres), and the error becomes non-uniform on a graded mesh. Should be
`cumsum(x.*dzCell) − 0.5*x.*dzCell`.

**Cost, and why it is study-only.** dt is a min over cells and dispersion is dz⁻², so
refining to √κ = 0.316 mm tightens the global bound **~556×**. A ~300-cell graded mesh is
~11× cheaper than uniform N = 3162 for the same interface resolution, but still **~1700×**
today's runtime — before the new dz⁻⁴ contribution that the κ fix introduces via `vbmax`.
Removing the dz⁻² restriction would need implicit dispersion (precedent: implicit osmosis,
`SSF.jl/src/simulate.jl:421-427`, added for exactly this reason) — explicitly out of scope.

**Bring forward if:** the κ-refinement convergence study shows resolved-interface results
differ materially from the sub-grid ones, i.e. if the interface actually matters to reported
quantities.

---

## Implicit Solver B

**Status:** investigated 2026-08-23, not started. This is where the timestep actually is —
Bailo is not (Phase 0: X = 1.000 at N=100/200/500).

### What Solver B is

Forward Euler on three blocks (`simulate.m:651-658`):

```
globalBiofilm += (dt/dz)(fluxIn - fluxOut)/porosity + dt*rhsBiofilm
globalFlowing += (dt/dz)(fluxIn - fluxOut)/porosity + dt*rhsFlowing
phiW          += (dt/dz)(fluxIn - fluxOut)/porosity + dt*rhsEnclosedWater
```

At N = 500: 1002 cells, 23 unknowns per cell (biofilm 13 = matrix 4 + enclosed particles 4 +
enclosed liquids 5; flowing 9; enclosed water 1), **23046 total**.

**Spatial coupling enters only through the fluxes**; reactions and exchange are entirely
cell-local. So the problem splits as follows -- but see the correction below on what the
dispersion row really depends on:

| part | structure | cost |
|---|---|---|
| dispersion (flowing only) | tridiagonal per component **only if phiFlowing is lagged** (see below) | 9 tridiagonal solves of size 1002 |
| advection (all blocks) | bidiagonal (upwind), components uncoupled | trivial |
| reactions + exchange | dense but LOCAL, 23x23 | 1002 independent small solves |

A fully coupled 23046 x 23046 implicit solve would be a serious build; the split makes it
routine. Nothing here needs a global nonlinear solve.

### CORRECTION (2026-08-23): flowing components DO couple through dispersion

An earlier version of this entry claimed "components do not couple through dispersion, so
this is one tridiagonal solve per component". That is wrong as stated. Dispersion acts on the
LOCAL concentration `g ./ phiFlowing`, and (`simulate.m:337-340`)

    phiMatrix   = sum(globalMatrix,2)/densityP;
    phiEnclosed = phiW + sum(globalEnclosedP,2)/densityP + sum(globalEnclosedL,2)/densityL;
    phiBiofilm  = phiMatrix + phiEnclosed;
    phiFlowing  = 1 - phiBiofilm;

so `phiFlowing` is one minus a **sum over all 13 biofilm-block components plus enclosed
water**. What is actually true:

- At FIXED `phiFlowing`, flowing component j does not see flowing component k. The 9
  independent tridiagonals are therefore valid **only because `phiFlowing` is lagged**.
- A genuinely implicit dispersion couples the flowing block to all 14 biofilm-block unknowns
  per cell, giving a block-banded system of 23 per cell -- not 9 scalar tridiagonals. The
  cheap structure is a property of the LAG, not of the operator.

**Why this matters beyond bookkeeping.** The sensitivity of the operator to the lagged
quantity is `d(1/phi_f)/d phi_f = -1/phi_f^2`. At `phi_f ~ 0.86` (supernatant, phi_b ~ 0.14)
that is ~1.35 and harmless. As the filter clogs, `phi_f -> 0` and it **diverges**, so the lag
degrades exactly in the regime the schmutzdecke work is about, and it compounds the
conditioning problem in blocker 7. This is a plausible mechanism for the predicted 14.5x
either failing to appear or appearing early in a run and breaking late.

Worth measuring directly: track `max(1/phi_f^2 * |d phi_b/dt| * dt)` over a run as a proxy for
the per-step error the lag introduces, and check whether it grows toward the end.

### What it is worth, from Phase 0 measurements (not estimates)

- **Dispersion is 88.6% of the binding region's sum at N = 500.** Removing it scales the two
  flowing regions by (1 - 0.886) and leaves regions 1-3 untouched, so the projected shares
  become matrix 0.0688, enclosed P 0.0693, **enclosed L 0.1345**, flowing P/L 0.1096 --
  i.e. **enclosed L binds next** and the gain is **7.43x** -> dt ~ 2.2e-05 d.
  (An earlier version of this entry said 14.5x, taken from the MATRIX margin 0.0688. Wrong
  region: the gain is set by the LARGEST surviving region, not the smallest.)
- **What is stiff after that: reactions.** `w_a` is identically zero in the enclosed-liquid
  region, so the survivors there are `w_v = 2*vbmax`, `w_b` (liquid transfer, plus osmosis
  `1/tau` unless ImplicitOsmosis), and `w_s = max|L_b| + max|L_e|` with
  `L = -sum(sigma_L .* mu .* X ./ (S + K_CFL))`. That last one is the Monod
  liquid-consumption bound and goes stiff as `S -> 0`, i.e. on substrate depletion --
  dissolved oxygen being the obvious candidate given the measured anoxia (32% of bed samples
  below 1 mg/L against a 7.0 mg/L day-mean). NOT YET CONFIRMED: the instrumentation records
  term shares only for the argmax region, so enclosed-L's composition is unmeasured. Job
  3531640 reports it directly.
- Then **advection binds**: dz/q = 2.78e-04 d = 24 s, i.e. **93x** current MaxDt in total.
- The 1/24 d target is a further **150x** beyond that and needs implicit advection, which
  means running 150x past the advective CFL. That is an ACCURACY problem, not a stability
  one: the transport would be badly smeared. **1/24 d is likely not reachable with a
  meaningful solution.**

So: implicit dispersion is the single biggest, cheapest win. Implicit reactions come next.
Implicit advection is last and probably not worth it.

### Blockers, concrete

1. **The Liebig minimum is non-smooth.** `evaluateReactions` (`simulate.m:796`) takes
   `min(mon_term, [], 2)` across Monod terms. Implicit reactions therefore need a
   **semismooth** Newton, freezing the argmin as an active set. The good news: this is the
   same machinery already written for the Bailo prototype
   (`analysis/prototypes/bailo/bailoCH1D.m`), which freezes active sets for `sign(u)` and for
   the clipping in `(1+x)^+`. Reusable.

2. **Reaction rates are products of state.** `rx = phi.*mu.*I.*monod.*product` with
   `product = local(:, orders)`, so an implicit reaction step is a genuine polynomial
   nonlinear system per cell -- fine, but the 23x23 Jacobian must be assembled per cell.

3. **Velocity coupling back to Solver A.** `velBiofilm` comes from Solver A, which needs
   `phi_b` from Solver B. Fully implicit means coupling them; lagging keeps the spatial part
   linear but the lag itself imposes a dt ceiling that **must be measured, not assumed** --
   this is the most likely way the 14.5x quietly fails to materialise.

4. **Implicitness does NOT buy positivity.** The negativity guard (`simulate.m:602,614`)
   would still be needed. Implicit Euler on Monod consumption is positivity-preserving under
   conditions, but the exchange terms (attachment, detachment, transfer) and the
   stoichiometric couplings are not obviously so. Same class of problem as the CH bounds, and
   NOT solved by the same trick: Bailo's upwind mobility argument applies to advective fluxes,
   not to reaction sources.

5. **State-dependent dispersion coefficient.** `dispersionStrength = |velFlowing|*(1-phi_b)`.
   Lagging it keeps dispersion linear and tridiagonal; that is almost certainly acceptable but
   is another lag whose dt ceiling is unmeasured.

6. **Upwind switches** `max(0,v)` / `min(0,v)` are constants if velocities are lagged, so they
   are only a problem for fully implicit advection -- another reason to leave that last.

7. **Clogging stiffness persists.** The `1/(1 - maxPhib)` factor blows up as phi_b -> 1.
   Implicit treatment removes it from the CFL but not from the operator's conditioning.

### Precedent in the tree

`ImplicitOsmosis` (`simulate.m:30, 448, 654`) already does exactly this for one term, as
`rhsEnclosedWaterB = rhsEnclosedWater./(1 + dt*kosm)` -- a scalar relaxation, the trivial case.
It is a precedent for the pattern and for gating it behind a default-false option.

### Suggested order

1. Implicit dispersion only, behind an option, velocities AND `phiFlowing` lagged. Linear,
   tridiagonal, 9 solves -- cheap only because of that lag, per the correction above.
   **Measure the achieved dt against the predicted 7.43x** -- if it falls short, the
   lag in (3)/(5) is why, and that is worth knowing before any further work.
2. Implicit reactions + exchange, local semismooth Newton reusing the Bailo Jacobian pattern.
3. Implicit advection only if 1/24 d is still wanted after seeing what (1) and (2) give, and
   with the accuracy cost quantified first.

### Sketch: dispersion without the phiFlowing lag

Two levels. **The cheap one is probably sufficient and should be tried first.**

#### Level A — use phi_f^{n+1}, which is already available. No lag, no extra cost.

The biofilm block carries **no dispersion at all** (`w_a = [0 0 0, ...]`, `simulate.m:562`;
`fluxBiofilm` is pure upwind advection). And in the present scheme every source and flux for
that block is evaluated at time n. Therefore

    globalBiofilm^{n+1}  and  phiW^{n+1}   are fully determined by time-n data,

so `phi_f^{n+1} = 1 - phiBiofilm^{n+1}` is **computable before the flowing solve**. Reordering
the step to

    1. advance globalBiofilm and phiW explicitly      -> phi_f^{n+1} known
    2. solve the flowing dispersion implicitly USING phi_f^{n+1}

removes the phi_f lag entirely, keeps 9 independent tridiagonals, and costs nothing. This is
not "fully implicit", but it eliminates the specific approximation that worried us.

What remains lagged at level A is `S_f = |v_f| * phi_f^face`, because `velFlowing` depends on
`velBiofilm` from Solver A. The `phi_f^face` half can also be taken at n+1; only `|v_f|`
genuinely needs Solver A.

#### Level B — genuinely implicit, needed only if the biofilm block also goes implicit

Write the per-cell state as `y` (14: 13 biofilm columns + enclosed water) and `g` (9 flowing).
`phi_f` enters as an **affine functional of y**:

    phi_f = 1 - a' y,     a = [1/rho_P (x8 particle cols), 1/rho_L (x5 liquid cols), 1 (water)]

The dispersive flux is

    F_face = |v_f| * phi_f^face * alpha_j * ( g_+/phi_f+ - g_-/phi_f- ) / dz

**Key structure: y enters ONLY through the scalar phi_f.** So the coupling of the dispersion
block to all 14 biofilm unknowns is **rank one per cell**. The clean formulation introduces
`phi_f` as one auxiliary unknown per cell with the linear constraint `phi_f + a'y = 1`,
reducing the dispersion system from 23 to **10 unknowns per cell** (9 flowing + 1).

Newton pieces:

    d(1/phi_f)/dy_k = + a_k / phi_f^2          (since d phi_f/dy_k = -a_k)
    dc_{i,j}/dy_k   = g_{i,j} * a_k / phi_f^2

so `dDISP/dy = (column vector depending on g) * a'` -- rank one, and cheaply handled by
Sherman-Morrison or by carrying the auxiliary unknown. The `dDISP/dg` block is exactly the
tridiagonal already implemented, with `1/phi_f` inside.

Cost: block-tridiagonal with 10x10 blocks over ~1002 cells per Newton iteration, versus 9
scalar tridiagonals at level A. Perhaps 5-10x level A per iteration, times 2-4 iterations.

#### Correction to the sensitivity estimate given earlier

The earlier note said the lag error scales as `-1/phi_f^2` and diverges as the filter clogs.
That overstates it, because `dispersionStrength` already carries a factor `phi_f`:

    dispersionStrength = abs(velFlowing) .* (1 - phiBiofilmBoundaries)      % == phi_f at faces

so with `phi_+ = phi + delta/2`, `phi_- = phi - delta/2`, `phi_face ~ phi`,

    F ~ |v_f| * alpha * [ (g_+ - g_-) - (g/phi_f) * delta ] / dz,     delta = grad(phi_f)*dz

The leading term is **phi_f-independent** -- the `phi_f^face` multiplying and the `1/phi_f`
dividing largely cancel. The residual sensitivity is `~ g * grad(phi_f) / phi_f`, i.e. first
order in `1/phi_f` and proportional to the GRADIENT of `phi_f`, not to `phi_f` itself. So the
lag is harmless wherever `phi_f` is smooth, and only bites where it is both small and steep --
the schmutzdecke front. Still worth avoiding, but level A avoids it for free, and the earlier
`1/phi_f^2` claim was too pessimistic.

### Measured: implicit dispersion at N=200, 3 d — and a corrected scaling argument

First production-regime measurement (job 3531640, task 1):

| | |
|---|---|
| steps | 302170 -> 106116, **2.85x reduction** |
| wall time | 234 s -> 125 s, **1.87x** (the 9 tridiagonal solves cost real time) |
| accuracy, explicit vs implicit at common dt | **1.56e-06** |
| binds now | **enclosed L**, 72% of steps |
| its composition | `w_v/dz` **0.323** | `w_a` 0.000 | `w_b` 0.018 | `w_s` **0.658** |

So reactions do become the stiff term, as projected. Note the gain (2.85x) came in BELOW the
4.97x projected from region shares -- the linear projection assumed time-mean term shares
compose exactly, which they do not. Treat the 7.43x quoted for N=500 as an upper bound.

**Corrected scaling argument.** An earlier version of this reasoning claimed that making
reactions implicit "trades N^1 for N^2 and loses" at fine mesh, because `w_s` is dz-free
(cost ~ N) while `w_v/dz` is not (cost ~ N^2). That is wrong. The N^1 regime is TRANSIENT:
`w_v/dz` grows without bound under refinement and `w_s` does not, so advection overtakes
reactions **on its own** at N ~ 407 (from A/S = 0.323/0.658 at N=200). Past that the explicit
scheme is also N^2, so there is no divergence to lose against.

| N | w_v share | w_s share | binds | implicit-reaction gain |
|---|---|---|---|---|
| 200 | 0.323 | 0.658 | reactions | 2.93x |
| 407 | 0.496 | 0.490 | crossover | 1.96x |
| 500 | 0.544 | 0.443 | advection | 1.80x |
| 800 | 0.656 | 0.334 | advection | 1.50x |
| 1615 | 0.794 | 0.200 | advection | 1.25x |

Implicit reactions are therefore **diminishing but never negative**. The argument against
them is not that they backfire, it is that 1.25x at N=1615 is a poor return on a per-cell
23x23 semismooth Newton over a non-smooth Liebig min, executed 1002 times per step.

**What the table actually points at:** at N >= 500 the binding term is ADVECTION, reaching
79% at N=1615. Implicit upwind advection is bidiagonal and structurally simpler than
reactions. The obstacle there is accuracy, not implementation -- `dz/q` at N=1615 is
8.6e-05 d, and exceeding it smears the transport.

### CLOSED: implicit Solver B is stage 1 only

Decision 2026-08-23, Jaime: *"explicit advection is fine, dt ~ dz is tolerable."*
(As spoken this was "dt ~ 1/dz", corrected immediately after; the advective CFL is
dt <= dz/q, so dt is PROPORTIONAL to dz. All the analysis below already used dt ~ dz.)

That settles the three-stage plan above:

| stage | verdict |
|---|---|
| 1. implicit dispersion | **BUILT.** Removes the dz^-2 term -- the intolerable one. 2.85x steps / 1.87x wall at N=200, accuracy 1.56e-06. |
| 2. implicit reactions | **NOT WORTH BUILDING.** Diminishing: 2.93x at N=200, 1.80x at N=500, 1.25x at N=1615, for 1002 per-cell 23x23 semismooth Newton solves per step over a non-smooth Liebig min. |
| 3. implicit advection | **NOT NEEDED.** dt ~ dz accepted. |

With dt ~ dz, cost scales as N^2, and implicit dispersion contributes a 1.87x wall gain at
fixed N:

| target | cost vs today's N=500 explicit |
|---|---|
| N=500, kappa=1e-6 (current) | 0.20x |
| N=800 | 0.51x |
| **N=1615, kappa=1e-7 RESOLVED (dz/ell = 1.00)** | **2.1x** |

(Updated once N=500 was measured. The first version of this table used the N=200 wall gain of
1.87x; the measured N=500 gain is **4.97x**, because the per-step tridiagonal overhead is
roughly fixed and so dilutes as the mesh grows. N=1615 is therefore ~2.1x a current N=500
run, not 5.6x.)

**Consequence worth acting on: the PUBLISHED kappa = 1e-7 is now affordable at a genuinely
resolved mesh.** Previously N=1615 cost ~34x an N=500 run (N^3 under explicit dispersion),
which is why kappa=1e-6 at N=500 was recommended instead. At 5.6x it is a routine cluster
job, so the manuscript need not defend a 10x change in kappa. Advective CFL there is
dz/q = 8.6e-05 d = 7 s.

Remaining open item, if dt ever becomes painful again: advection is 79% of the binding sum
by N=1615, and implicit upwind advection is bidiagonal -- structurally simpler than
reactions. The obstacle is accuracy (exceeding dz/q smears transport), not implementation.

### Measured: implicit dispersion at N=500, 3 d — 7.11x, and the crossover confirmed

Job 3531640 task 2, production regime:

| | N=200 | **N=500** |
|---|---|---|
| steps | 302170 -> 106116, 2.85x | 1448486 -> 203600, **7.11x** |
| wall | 234 -> 125 s, 1.87x | 2289 -> 460 s, **4.97x** |
| accuracy vs explicit, common dt | 1.56e-06 | 6.10e-05 |
| binds now | enclosed L, 72% | **enclosed L, 92%** |
| `w_v/dz` : `w_b` : `w_s` | 0.32 : 0.02 : 0.66 | **0.53 : 0.01 : 0.46** |

**7.11x against the 7.43x projected** -- within 4%. The projection method had looked unreliable
after N=200 came in at 2.85x against a projected 4.97x, but it holds at the mesh that matters.

**The advection/reaction crossover is confirmed and lands where predicted.** From the N=200
shares the crossover was projected at N ~ 407; at N=500 `w_v/dz` = 0.53 now exceeds
`w_s` = 0.46. A forward prediction from independent data, which validates the scaling table
used for the kappa=1e-7 costing.

Consequence: **reactions are NOT the stiffest term at production resolution.** They were at
N=200; at N=500 advection is. Implicit reactions would buy ~1.8x here and less under
refinement -- confirming they are not worth building.

Accuracy caveat unchanged: 6.10e-05 is 40x the N=200 figure, and this run PREDATES the level-A
reordering, so part of it may be the phi_f lag rather than the splitting. Separate the two
before relying on this at N=1615.
