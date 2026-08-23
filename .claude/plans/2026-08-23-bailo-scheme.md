# Bailo2023 bound-preserving Cahn–Hilliard scheme — implementation plan

**Paper.** Bailo, Carrillo, Kalliadasis & Perez, *Unconditional Bound-Preserving and
Energy-Dissipating Finite-Volume Schemes for the Cahn–Hilliard Equation*
(arXiv:2105.05351v2), `documents/Bailo2023.pdf`. §2 is the 1-D scheme, Eq. (2.1a)–(2.1i).

**Goal.** Decide, on evidence, whether to replace Solver A with Bailo's scheme — and if yes,
have a validated standalone prototype ready to integrate. "Done" for this plan is a decision
plus a prototype, **not** an integration.

---

## What was established before writing this plan

Three facts checked against the paper and the code, not assumed.

**1. The transform is exact, and it is not just a relabelling of the mobility.** Bailo's
phase parameter lives in (−1,1); our volume fraction lives in (0,1). The linear map is
φ_B = 2u − 1, u = (1 + φ_B)/2. Pushing it through *every* term — not just the mobility —
gives three factors that must be tracked together:

    dφ_B/dt = 2 du/dt          (1 + φ_B)(1 − φ_B) = 4u(1 − u)          ∂²φ_B/∂z² = 2 ∂²u/∂z²

and, by the chain rule on the energy variation, μ = δF/δu = 2 δF/δφ_B = 2ξ, so ξ_z = μ_z/2.
Imposing that normalisation and matching term by term against our own system
`u_t = ∂_z[ζ₀u(1−u) ∂_z μ]`, `μ = Ψ′(u) − κ ∂²u/∂z²`:

| Bailo | ours |
|---|---|
| M₀ | **ζ₀** |
| ε² | **κ/4** |
| H(φ_B) | **Ψ((1 + φ_B)/2)**, i.e. Ψ′(u) = 2H′(2u − 1) |

**Verified numerically**, both formulations assembled on a 2000-cell grid: max relative
difference 1.3e-10 (finite-difference truncation of the nested derivatives). The obvious
naive readings are badly wrong — M₀ = ζ₀/4 with ε² = κ gives **75% error**, and ε² = κ with
M₀ = ζ₀ gives 1.8e-04. The factor of 4 in the mobility does *not* survive as a factor of 4 in
M₀: it is cancelled by the two factors of 2 from dφ_B/dt and ξ_z.

**Implement natively in (0,1); do not round-trip through φ_B.** Transforming our data into
Bailo's variable, running their scheme, and transforming back invites exactly the factor
errors above at every interface, and makes the prototype read in a variable nobody else in
the model uses. The scheme rewrites cleanly. In particular the two-point upwind mobility
becomes

    M(x, y) = ζ₀ · (x)⁺ · (1 − y)⁺

(checked against M₀(1+x_B)⁺(1−y_B)⁺/4 at interior, out-of-range and negative arguments: exact
agreement), and the bound-preservation argument transfers verbatim with 0 and 1 as the
bounds — if u_k > 1 then (1 − u_k)⁺ = 0, and if u_j < 0 then (u_j)⁺ = 0, killing the same two
flux terms in the same contradiction proof.

The velocity also maps onto something we already compute. Bailo's u_{i+1/2} = −(ξ_{i+1} − ξ_i)/Δx
is our v_b up to the mobility placement: `run_biofilm.m` forms
`v_b = q − ζ₀(1−u)(μ_{i+1} − μ_i)/Δz`, so the flux u·v_b = ζ₀u(1−u)(−∇μ) is already the
degenerate form Bailo requires. Only the *upwinding* of that mobility is new.

**2. Bound preservation does NOT depend on the potential.** This is the load-bearing fact.
The contradiction proof (§2.1, p. 6) sums the scheme over a contiguous group of cells with
φ > 1 and kills the two dangerous flux terms using only `(1 − y)⁺ = 0`, a property of the
*upwind mobility form* M(x,y) = M₀(1+x)⁺(1−y)⁺. The double-well structure of H is never
used. Our potential is therefore no obstacle to the main prize.

**3. Our potential is nearly convex, and is not a double well.** With
Ψ′(φ) = φ²(φ − 1.5ζ₁), ζ₁ = 0.01:

- Ψ″(φ) = 3φ(φ − ζ₁), minimum **−7.5e-05** at φ = ζ₁/2 = 0.005; spinodal band only
  0 < φ < 0.01.
- So the convex-splitting constant needed is **a = 7.5e-05**: H_c = Ψ + aφ²/2, H_e = aφ²/2
  are both convex on [0,1]. Trivially satisfiable.
- Ψ′ has roots at φ = 0 (double) and φ = 1.5ζ₁ = 0.015. Ψ(0.015) = −4.2e-09, Ψ(1) = 0.245.
  There is **one shallow interior minimum and no second well near φ = 1.**

Consequence, stated plainly: **the energy-dissipation half of Bailo buys us almost nothing**
— our potential is convex to within 7.5e-05, so a convex splitting is nearly vacuous. The
reason to do this is bound preservation, and bound preservation alone. Any write-up should
say so rather than claiming a discrete energy law as a headline benefit.

---

## Two obstacles the paper does not cover

Both are real and neither is mentioned in the parked `PLANS.md` entry.

**A. We have a reaction source; Bailo does not.** `run_biofilm.m:349` solves with
`RHS = [un + dt*Rbio(1:n0); DPSI]`. Bailo's Eq. (2.1a) is a pure conservation law,
∂φ/∂t + ∇·F = 0, and *both* proofs (mass conservation and boundedness) rely on that. Biofilm
growth can push φ past 1 no matter how well-behaved the flux is. Bound preservation for the
combined system is therefore **not** inherited.

Resolution to prototype: Strang/Lie operator splitting — a bound-preserving CH substep, then a
reaction substep with its own positivity handling. The reaction substep needs its own bound
argument (or a cap that is a *physical* clogging event rather than a numerical abort).

**B. We have convection; Bailo does not.** Our CH block carries the `S_q` convection matrix
(inflow q̂ through the supernatant). Adding advection breaks the boundedness proof unless the
advective flux is itself upwinded and folded into the same contradiction argument. The scheme
already upwinds the cohesive flux, so extending the same treatment is plausible — but it must
be *proved or numerically demonstrated*, not assumed.

These two determine whether Bailo's "unconditional" claim survives contact with our model. If
they don't, the honest outcome is "bound-preserving CH substep inside a scheme that is not
unconditionally bound-preserving overall" — still useful, but a weaker claim.

---

## Phase 0 — Is it worth it? (measurement first, ~half a day)

Do not write scheme code until this is answered. Instrumentation only; no change to the
numerics of `simulate.m`.

**What the CFL bound actually looks like** (`simulate.m:491`):

    dt_CFL = cflFactor / max( w_v/dz + w_a/dz^2 + w_b + w_s )

Each `w_*` is a 5-vector over REGIONS — `[matrix, enclosed P, enclosed L, flowing P,
flowing L]` — summed within a region, then maxed across regions. So "which term binds" is a
share-of-sum inside the argmax region, not a max over terms. Two facts from the code:

- `w_v = [2*vbmax*[1 1 1], 2*vfmax*[1 1]]`, so **κ reaches dt through exactly one quantity,
  `vbmax = max|velBiofilm|`**, in regions 1-3.
- `w_a = [0 0 0, ...]` — dispersion is identically zero in the biofilm/enclosed regions, so if
  the matrix region binds, the competition is only `w_v/dz` vs `w_b` vs `w_s`.

**The ceiling on any possible speedup.** `v_b = q − ζ₀(1−φ)∇μ`. Bailo treats the cohesive term
implicitly but cannot remove `q`. The best it can do is delete the cohesive share of `vbmax`.
If `vbmax` is mostly advective, Bailo buys no speed at all, however good the scheme is.

**Four measurements, in order; each can terminate the investigation.**

0a. **Is dt even CFL-bound?** Fraction of steps with `dt_CFL < adaptiveMaxDt`
    (`manuscriptExperiments` passes `MaxDt=3e-6`). If dt sits on the cap, every downstream
    number is meaningless. This exact trap produced a meaningless κ-fix dt-cost measurement
    on 2026-08-23 — do not repeat it.
0b. **Which region, which term?** For the argmax region, the fractional contribution of
    `w_v/dz`, `w_a/dz^2`, `w_b`, `w_s`, as a histogram over steps.
0c. **How much of `vbmax` is cohesive?** Recompute `vbmax` at the same state with the κ term
    zeroed. The ratio is the counterfactual: *if the CH restriction vanished entirely and for
    free, dt would rise by X×.* **This is the decision number.**
0d. **Does the negativity guard trip?** Count `Flag = "BIOFILM"` aborts
    (`simulate.m:555-578`) across production configurations.

**Regime requirement.** Run from the mature state of the 20 d run (job 3531519), NOT a clean
filter. On a clean start φ_b ≈ 0, the cohesive term is inert and 0c returns X ≈ 1 for reasons
that have nothing to do with the scheme.

**Decision gate.**
- 0c gives X ≫ 1 → the speed case is real; continue to Phase 1.
- 0c gives X ≈ 1 but 0d shows the guard tripping → continue, but for **bound preservation
  only**. Say so explicitly, and drop any speed claim from the write-up.
- 0c gives X ≈ 1 and 0d shows no aborts → re-park, with the numbers recorded.

Deliverable: `analysis/probes/probeCflBudget.m` and a short note in `.claude/decisions/`.

## Phase 1 — Standalone prototype, the paper's own test case

Scratch MATLAB in `analysis/prototypes/bailo/` (not `src/`), no reaction, no convection.
Implement (2.1a)–(2.1i) **verbatim in Bailo's own φ ∈ [−1,1]** for this phase only — the
validation case and its exact steady state are stated in that variable, and reproducing them
literally is what makes Phase 1 an independent check of our implementation. Phase 2 then
rewrites the same code natively in (0,1). Keep the two in one file behind a variable-convention
switch so the Phase 1 gates can be re-run against the native form as a regression.

- upwind flux (2.1b) with M(x,y) = M₀(1+x)⁺(1−y)⁺
- velocity (2.1c), semi-implicit ξ (2.1e) with H_c implicit / H_e explicit and Δφ implicit
- no-flux via F_{1/2} = F_{M+1/2} = 0 and ghost values (2.1i)
- wetting term W set to zero (we have no contact-angle physics; the boundary is z = 0 Neumann)

Nonlinear solve: Newton on the M-cell system for φⁿ⁺¹, analytic Jacobian, tol 1e-12. The
paper does not specify its solver, so this is our choice and must be reported.

**Validation — §4.1, which has an exact answer.** Deep-quench logarithmic potential
H(φ) = (1 − φ²)/2, degenerate mobility, datum φ₀ = cos(x/ε) − 1 on |x| ≤ πε/2 else −1, domain
[−3πε/2, 3πε/2], run to T = 20ε². Exact steady state (4.1):

    φ∞(x) = (1/π)(1 + cos(x/ε)) − 1 on |x| ≤ πε, else −1

Run at ε = 1, 0.1, 0.01, 0.001 — the paper insists on the range, and so should we.

**Gates (all must pass):**
| gate | criterion |
|---|---|
| mass | Σφⁿ conserved to machine precision |
| bounds | \|φ\| ≤ 1 at every cell, every step, **including at large Δt** |
| energy | F_Δ[φⁿ⁺¹] ≤ F_Δ[φⁿ] monotonically, per (2.2)/(2.3) |
| accuracy | second order in Δx against φ∞, at every ε |
| unconditionality | bounds and energy hold as Δt is raised by ≥ 100× past the explicit limit |

The last gate is the whole point; a prototype that only works at small Δt has proved nothing.

---

## Phase 1 results (2026-08-23) — prototype built, three gates pass, two open

`analysis/prototypes/bailo/bailoCH1D.m` implements (2.1a)-(2.1i) in Bailo's own
phi in [-1,1]; `validateBailo.m` runs the gates. Deep-quench case of section 4.1,
H = (1-phi^2)/2, so Hc = 0 and He = phi^2/2 - 1/2.

**Nonlinear solve.** The paper specifies none. Used here: Picard with the mobility AND
the upwind direction lagged, making each iterate a linear pentadiagonal solve. Lagging is
legitimate at convergence -- the converged iterate satisfies the fully implicit equations,
which is what the proofs need.

**Gates 1-3 PASS**, across the paper's range of eps:

| eps | mass drift | bound violation | max energy increase | Picard iters |
|---|---|---|---|---|
| 1    | 8.62e-13 | 3.95e-11 | 1.30e-12 | 33.6 |
| 0.1  | 1.45e-12 | 3.55e-11 | 2.41e-13 | 28.9 |
| 0.01 | 2.69e-12 | 4.69e-11 | 1.80e-14 | 31.6 |

Bound violations sit at the same order as the mass drift and the operator's conditioning,
so they read as solve roundoff, not structural. Energy is monotone to roundoff.

**Gate 4 (second order) NOT YET ESTABLISHED.** Convergence to the exact steady state (4.1):

| M | dx | err | order |
|---|---|---|---|
| 50  | 1.885e-02 | 7.7435e-04 | - |
| 100 | 9.425e-03 | 2.0710e-04 | 1.903 |
| 200 | 4.712e-03 | 5.0453e-05 | 2.037 |
| 400 | 2.356e-03 | 1.3393e-05 | 1.913 |
| 800 | 1.178e-03 | 5.7400e-06 | **1.222** |

Second order over three refinements, then it degrades. Most likely the error floor: T is
fixed at 20*eps^2 and dt at 0.1*eps^2, so `err vs phi_inf` mixes spatial discretisation with
leftover transient and first-order temporal error. Refining the mesh shrinks the first only.
Test running: hold M = 800 and grow T. **Not yet confirmed.**

**Gate 5 (unconditionality) NOT ESTABLISHED -- and this is the important one.**

| dt | mass | bounds | max dE | Picard | err |
|---|---|---|---|---|---|
| 0.001 (1x)    | 1.45e-12 | 3.55e-11 | 2.41e-13 | 28.9 | 5.045e-05 |
| 0.01 (10x)    | 1.49e-11 | 5.18e-09 | 0 | **2000 (cap)** | 2.145e-04 |
| 0.1 (100x)    | 8.91e-12 | 1.05e-10 | 0 | **2000 (cap)** | 7.370e-02 |
| 1 (1000x)     | 3.67e-10 | 0 | 0 | **2000 (cap)** | 1.536e-01 |

Bounds and energy do hold at 1000x dt -- but **Picard hit its iteration cap at every dt
beyond 1x**, so those runs are not solving the implicit system at all, and the zeros in the
energy column are consistent with a stalled iterate rather than a dissipating one. The
unconditional-stability claim is therefore untested. The weak link is the ITERATION, not the
scheme.

**Consequences for the plan.**
- Phase 3's cost question is already partly answered and the answer is bad: ~29-34 linear
  solves per step against the present scheme's ONE. Even at 1x dt this is ~30x the Solver A
  cost, and Phase 0 established there is no dt to be won back (X = 1.000).
- Picard must be replaced by Newton (or Anderson-accelerated / damped iteration) before
  Gate 5 can be attempted. The upwind max/min makes the residual non-smooth, so this wants
  a semismooth Newton rather than a plain one.

## Phase 2 — Substitute our potential and mobility (native (0,1) form)

Same prototype, still standalone, now rewritten in u ∈ [0,1] — no change of variable at
runtime. The scheme, stated natively:

    (u_i^{n+1} − u_i^n)/Δt + (F_{i+1/2}^{n+1} − F_{i−1/2}^{n+1})/Δz = 0                (2.1a)

    F_{i+1/2} = M(u_i, u_{i+1})(v_{i+1/2})⁺ + M(u_{i+1}, u_i)(v_{i+1/2})⁻              (2.1b)
    M(x, y)   = ζ₀ (x)⁺ (1 − y)⁺                                                       (2.1d)
    v_{i+1/2} = −(μ_{i+1} − μ_i)/Δz                                                    (2.1c)

    μ_i^{n+1} = Ψ_c′(u_i^{n+1}) − Ψ_e′(u_i^n) − κ (Δ_h u)^{n+1}_i                      (2.1e)
    Ψ_c = Ψ + a u²/2,   Ψ_e = a u²/2,   a = 7.5e-05                                    (§3 above)

all quantities at n+1 except Ψ_e′, with no-flux F_{1/2} = F_{M+1/2} = 0 and ghost values
u_0 = u_1, u_{M+1} = u_M. κ enters directly — **there is no ε² and no factor of 4 anywhere in
the native form**; the κ/4 appears only if one insists on writing it in Bailo's variable.

Re-run the Phase 1 gates. Expect bounds and mass to hold (they are potential-independent) and
energy to stay monotone (the split is valid, a = 7.5e-05). **Do not expect phase separation**
— with one shallow minimum at u = 0.015 and no second well, the solution relaxes toward 0.015
with the u = 1 bound enforced by mobility degeneracy alone. If phase separation appears,
something is wrong.

Then compare against the current scheme on identical initial data at κ = 1e-6, N = 500: same
steady state? same interface width ℓ = √(κ/Ψ″) ≈ 1.96 mm?

**Cross-check to run once, cheaply:** assemble the native RHS and the transformed-Bailo RHS on
the same grid and confirm they agree to truncation, as was done when writing this plan. That
pins the constants and catches a transform error before it contaminates the gates.

## Phase 3 — Cost

- Newton iterations per step to 1e-10, and cost per step vs the current *single* linear solve
  of the 2n₀×2n₀ block system.
- Largest Δt at which each scheme is stable **and** accurate (stability is not enough — an
  unconditionally stable scheme can still be first-order-in-time inaccurate).
- Net: cost per simulated day. This is the number that decides integration.

---

## The payoff that makes boundedness structural, not just diagnostic

Raised 2026-08-23. This reframes what bound preservation is *for*, and is the strongest
argument in the plan.

Today Solver A solves a 2n₀×2n₀ system and **throws away half the solution**: uⁿ⁺¹ is
discarded (`run_biofilm.m`, the `%un1 = sol(1:n0)` line is commented out) and only v_b is
kept. Meanwhile φ_b is reconstructed *diagnostically* downstream as a sum over biofilm
components (`manuscriptExperiments.m:425-435`):

    φ_b = W_enc/ρ_L + Σ_particles (X_matrix + X_enc)/ρ_P + Σ_liquids S_enc/ρ_L

Nothing enforces φ_b ≤ 1. The sum is emergent, and when it misbehaves the code *aborts*
(`simulate.m:555-578`).

**If Bailo makes φ_bⁿ⁺¹ ∈ (0,1) trustworthy, we stop discarding it and spend it.** The
constraint Σ_{i∈b} φ_i = φ_b then eliminates one component from the transport system, which is
reconstructed algebraically instead of transported. Two consequences:

- the sum constraint holds **exactly by construction** rather than emergently, and
- φ_b ∈ (0,1) becomes **structural**, so the negativity abort turns into a physical clogging
  event.

**Which component to eliminate.** The reconstructed component is a residual,
φ_elim = φ_b − Σ(others), so it absorbs every error in φ_b and in all other components. Its
*relative* error is therefore minimised by eliminating the **largest** component. Measured on
`reference_pathogen`:

| component | max value | role |
|---|---|---|
| **Water / enclosed** | **1.01e-01** | largest; physically slaved (β = 0.99) |
| HET / matrix | 5.02e-02 | |
| PHO / matrix | 5.02e-02 | |
| POM / matrix | 2.01e-02 | |
| PAT / matrix | 1.00e-02 | smallest, and the scientific quantity of interest |

**Eliminate enclosed water.** It is the largest, and `BiofilmPorosity = 0.99` already ties
water content to biomass, so it is the component closest to being slaved already — the
elimination is arguably removing a redundancy rather than imposing one.

**Never eliminate PAT.** It is the smallest component and the pathogen result the manuscript
is about; reconstructing it as a small difference of larger numbers is cancellation error
aimed directly at the headline output.

**Two conditions that must be verified before this is sound.**

1. **Consistency of the two routes to φ_b.** The elimination is only legitimate if Solver A's
   φ_b equation *is* the sum of Solver B's component equations. That needs (i) every biofilm
   component advected with the same v_b, so Σφ_i v_b = φ_b v_b, and (ii) Σ_i R_i over biofilm
   components exactly equal to the `Rbio` source used in the CH step
   (`run_biofilm.m:349`, `RHS = [un + dt*Rbio(1:n0); DPSI]`). Any mismatch is dumped **entirely**
   into the reconstructed component. This is currently untested — and untestable today,
   because one of the two routes is discarded. **Check it first, before Bailo:** run the
   existing CH φ_b solve alongside the diagnostic sum and measure their divergence over a run.
   If they already disagree materially, that is a defect worth knowing about independently.
2. **The residual can still go negative.** φ_b ∈ (0,1) with the other components ≥ 0 does *not*
   imply φ_b − Σ(others) ≥ 0. Eliminating the largest component makes it unlikely, not
   impossible; a sign monitor is required. Note that if *enclosed water* is the one that goes
   negative, that is physically meaningful — biofilm solids exceeding the available volume —
   and is a legitimate clogging signal rather than a numerical failure.

**Effect on the Phase 0 gate.** This weakens the dependence on measurement 0c. If the
architectural win is the objective, the speed counterfactual matters less and the case rests
on 0d plus condition 1 above. Record which justification is being invoked.

**Note this is not achievable without bound preservation.** One could reconstruct a component
from the current CH solve today, but φ_b would carry no guarantee and the residual could be
arbitrary. Bailo and this elimination are complementary: the first makes the second safe.

## Phase 4 — Integration decision (decision only; integration is a separate plan)

Write the outcome to `.claude/decisions/`. If integrating:

- Bailo replaces Solver A. Our architecture already discards uⁿ⁺¹ and keeps only v_b, which
  maps onto Bailo's u_{i+1/2} — so the interface to Solver B is unchanged in shape.
- Obstacles A and B above must have a demonstrated answer.
- **Superseded objection.** An earlier draft said bounding φ_b is "necessary but not
  sufficient" because it does not bound each species. That framed boundedness as merely
  diagnostic. The section above uses it instead — see *The payoff that makes boundedness
  structural*. What survives of the objection is narrower: the single reconstructed component
  needs a sign monitor.
- Domain note in our favour: the CH system covers cells 1..n₀ (supernatant) with homogeneous
  Neumann at z = 0, which matches Bailo's no-flux (2.1h)/(2.1i) exactly.

---

## Reuse — read before writing anything

`code/2d/ssf-ch-stokes-fenics` already implements this family: Eyre convex splitting
(`main.py`, `dpsi_ey = PSI_0*(dpsi_i(u) + dpsi_e(u_0))`), mobility upwinding around
`u_m = ρ_P/(2+γ)` (`main.py:190-240`), mass lumping, and `tests/test_potentials.py` enforcing
that the split sums back to `f`. Its dolfinx port reports monotone energy decay and
machine-precision mass conservation. That is a working reference for the splitting and the
upwinding, in-house.

## Files

| Path | Change |
|---|---|
| `analysis/prototypes/bailo/` | **new** — standalone prototype, Phases 1–2 |
| `analysis/probes/probeCflBudget.m` | **new** — Phase 0 CFL-term attribution |
| `src/@State/simulate.m` | Phase 0 instrumentation only; no scheme change in this plan |
| `.claude/decisions/` | Phase 4 outcome |
| `PLANS.md` | unpark the Bailo entry, point it here |
