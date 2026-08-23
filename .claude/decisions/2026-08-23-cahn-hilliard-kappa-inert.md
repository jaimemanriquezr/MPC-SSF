# The Cahn–Hilliard interfacial term was inert; κ had no effect on any output

**Date:** 2026-08-23
**Status:** fixed, operator gate passed

## The defect

`getCahnHilliardMatrices.m` assembled all three operators with `Rows ≤ n0` — `convection`
(`:27`), `diffusion` (`:34`) and `diffusionMobility` (`:41`). The `2n0 × 2n0` mixed system was
therefore block **upper**-triangular:

| block | was | should be |
|---|---|---|
| (1,1) | `I − DD − Δt·S` | `I − Δt·S` |
| (1,2) | `−Δt·D_λ` | `−Δt·D_λ` ✓ |
| (2,1) | **`0`** | **`−DD`** |
| (2,2) | `I` | `I` ✓ |

κ enters only through `DD`, which `simulate.m:101` placed in block (1,1) via `CH0 = Id - DD`.
With block (2,1) empty and block (2,2) a bare identity, the second block row read
`I·μ = Ψ'(uⁿ)`, so **μ = Ψ'(u) exactly** and κ never reached the chemical potential. The φ
half of the solve, `xCH(1:n0)`, is discarded at `:431` — which is intended
(`reports/numerical-method.md:62`: Solver A exists only to produce `v_b`) — so the κ term in
block (1,1) affected nothing either.

Net effect: **the interfacial/surface-tension term of the Cahn–Hilliard model was switched
off** for every result produced to date.

**The defect is inherited, not introduced by the port.** The legacy implementation has it
too: `slow-sand-filtration/src/@SDfilter/get_ch_matrices.m:38` clamps `D.rows` to `[1, n0]`
while `D_m.cols` gets `+ n0`, and `run_biofilm.m:161` forms `G_base = Id - DD` with
`RHS = [un + dt*Rbio; DPSI]`. So all three implementations — legacy MATLAB, MPC-SSF, and
SSF.jl — collapse to `mu = Psi'(u)`. It originates in the thesis/legacy code and was
faithfully ported twice. **Manuscript Fig. 7 and every published figure were produced with
κ inert.**

## Evidence

*Operator level, independent of any simulation state.* Assembling `lhsCH` exactly as
`simulate.m` does, with random mobility and random `u`:

```
block (2,1)  nnz = 0        block (2,2) == I exactly
DD carries kappa, nnz = 61, all rows <= n0
max| muCH - psi'(u) | = 0.000e+00
```

*Corroboration.* `reports/log-oat-report-2026-08.typ` ranks κ 27/27 with **exactly** zero
sensitivity — `I_rms = I_max = D_min = I_min = A_i = 0` — in all three scenarios. That was
read at the time as κ being sub-grid, but sub-grid resolution gives small-but-nonzero.
Exactly zero is the signature of an inert parameter.

(An earlier draft also offered "the goldens run at κ = 1e-2 and pass, which they could only
do while κ did nothing". That inference is **wrong** and is retracted: they still pass after
the fix, because their near-uniform initial state makes the κ term inactive either way. See
the goldens note under Verification.)

*A discarded test, recorded so it is not repeated.* The first attempt swept κ over nine
orders of magnitude on a 0.05 d run from a clean start and reported bit-identical output.
That proves nothing: at 0.05 d there is no biofilm in the supernatant, so φ_b ≈ 0, Δφ ≈ 0 and
λ = ζ₀φ(1−φ) ≈ 0 — the term would have been inert *even if correctly wired*. Jaime caught
this. Any test of a mechanism must run in a regime where the mechanism is active; where a
state-independent check exists, prefer it.

## The fix

Eq. (29) is `μ_j = −κ(Du_{j+1/2} − Du_{j−1/2})/Δz + Ψ'(u_j)`, i.e. `μ = −κ·Δ_h u + Ψ'(u)`.
`getCahnHilliardMatrices.m:32-37` already assembles exactly `DD = κ·L` with `L = −Δ_h`
(scaled by `Δz²/κ` it is the tridiagonal `[1 −1; −1 2 −1; …; −1 1]`, homogeneous Neumann at
both ends). The μ row is therefore `μ − κ·L·u = Ψ'(u)`, so `DD` belongs in block (2,1).

`diffusionMobility` already reaches block (1,2) by adding `+ n0` to its **Columns**, so the
convention lives in this file. The whole fix is `+ n0` on `diffusion.Rows`; `CH0 = Id - DD`
at `simulate.m:101` then yields `−DD` in block (2,1) and `I` on the diagonal, with no change
to `simulate.m` at all.

**Not changed: the z = 0 boundary condition.** `numerical-method.md:250-262` shows a
`−κ·u₁^{n+1}/Δz²` term in the μ right-hand side. That belongs to the Diehl2025 **Dirichlet**
condition, which manuscript Appendix A.1 records as deliberately superseded: "we have
instead used the Neumann boundary condition ∂φ_b/∂z = 0 at z = 0", because the Dirichlet form
"would need very small time steps". The index clamping in `diffusion.Rows` already implements
that Neumann condition. `numerical-method.md` is transcribing the older scheme and should be
annotated.

## The legacy implementation carries the identical defect, and is now fixed too

The defect is **inherited**, not introduced by the port. `slow-sand-filtration` is the
published MATLAB original that MPC-SSF and SSF.jl descend from, and
`src/@SDfilter/get_ch_matrices.m:38` had the same clamp:

```matlab
D.rows   = min(max(1, dd_rows - 1), n0);          % block (1,1) -- kappa discarded
D_m.cols = cols(:) + (0:n0-2) + n0;               % block (1,2) -- correct
```

`run_biofilm.m:155-161, 347` assembles and solves it exactly as `simulate.m` does
(`G_base = Id - DD`, `G = G_base - dt*(S_q + D_lambda)`, `RHS = [un + dt*Rbio; DPSI]`,
`sol = G\RHS`), with `un1 = sol(1:n0)` commented out and only `mu = sol(n0+1:end)` used to
form `v_b`. Same structure, same consequence.

**Measured before changing anything** (published `SDparameters.mat` model, `zeta_0 = 1e2` as
`srun_biofilm.m` sets it, `add_cells(20)` → n0 = 21, κ = 1e-7, random `u` and random
mobility — a state-independent operator check, not a simulation):

| placement | block (2,1) nnz | block (2,2) == I | `max\|μ − Ψ'(uⁿ)\|` |
|---|---|---|---|
| as-is (rows clamped to 1..n0) | **0** | yes | **0.000000e+00** |
| with `+ n0` | 61 | yes | 5.886212e-06 |

and on the only quantity Solver A produces,
`max|v_b(fixed) − v_b(as-is)| = 1.458859e-02` (relative 6.516e-04). So the fix is
demonstrably needed here as well, and the check was run *before* the edit, per instruction.

Applied as a one-line change (`+ n0` on `D.rows`) plus a comment block, with no adjacent
cleanup and no style change — the workspace `CLAUDE.md` keeps this repo legacy and
snake_case, permitting changes only where they let it reproduce the newer implementations'
runs, which is exactly what parity of the CH scheme is for. Post-fix: `checkcode` clean, and
`run_biofilm` at `add_cells(30)`, `time_run = 0.25` returns `flag OK` with all values finite
(the single `Inf` in `flowing/WATER` is a pre-existing 1×1 sentinel, not a solver value).

**Consequence.** Every published figure — manuscript Fig. 7 included — was produced with
κ inert, in the legacy code as well as in MPC-SSF. The Appendix A.2 criterion Δz < √κ was
never violated at N = 500; it was *moot*.

## Also fixed, since both lied rather than failed

- `SandFilter.m` `GridSize` returned `mean(diff(Boundaries))`, which telescopes to
  `(last−first)/N` — silently wrong on any non-uniform mesh, with ~60 call sites carrying on
  regardless. It now errors on a non-uniform mesh; `CellWidths` and `CenterSpacings` added
  for the graded-mesh work parked in `PLANS.md`.
- `GridZero` used `find(...)` with **no `,1`**, returning a *vector* whenever more than one
  centre fell within half a mean width of zero — as happens on any mesh refined near z = 0 —
  and that vector then sized the CH system. Now `find(Boundaries >= 0, 1) - 1`, with an
  error if no such cell exists. The duplicate of both in
  `getCahnHilliardMatrices.m:13-14` is removed in favour of the filter's accessors.

## Verification

**Gate 1 (operator, state-independent) — passed.** Block (2,1) `nnz` 0 → 61; block (2,2)
still `I`; `max|muCH − Ψ'(u)|` 0 → 4.58 at κ = 1e-2. Scaling in κ:

| κ | max\|μ − Ψ'\| | ratio |
|---|---|---|
| 1e-8 | 4.628995e-06 | — |
| **1e-7** | 4.628994e-05 | **9.999999** |
| 1e-6 | 4.628990e-04 | 9.99999 |
| 1e-5 | 4.628949e-03 | 9.99991 |
| 1e-4 | 4.628539e-02 | 9.99911 |
| 1e-3 | 4.624442e-01 | 9.99115 |

Linear in κ as κ → 0, with deviation growing smoothly — as expected, since
`μ = Ψ'(uⁿ) − κ·L·u^{n+1}` and `u^{n+1}` is itself weakly κ-dependent through the coupled
solve. An earlier draft of the plan demanded exact linearity; that criterion was too strong.

**Gate 2 (manufactured solution) — passed, on cosmos (job 3531490).** A smooth
`u(z) = 0.45 + 0.30·cos(3πζ)` on the CH domain, chosen so `u' = 0` at both faces — otherwise
the homogeneous Neumann rows carry an O(1) consistency error unrelated to the operator being
tested. The assembled block (2,1) is applied to `u` and compared against the exact `−κ·u''`:

| N | dz | err interior | order | err all | order |
|---|---|---|---|---|---|
| 20 | 4.878e-02 | 4.128e-08 | — | 4.128e-08 | — |
| 80 | 1.242e-02 | 2.963e-09 | 1.964 | 2.963e-09 | 1.964 |
| 320 | 3.120e-03 | 1.908e-10 | 1.991 | 1.908e-10 | 1.991 |
| 1280 | 7.809e-04 | 1.201e-11 | 1.998 | 1.201e-11 | 1.998 |
| 2560 | 3.906e-04 | 3.006e-12 | **1.999** | 3.006e-12 | **1.999** |

Clean second order over seven refinements at κ = 1e-7, and again at κ = 1e-2 (order 1.995
over six). The interior and full-domain norms are **identical at every N**, so the index
clamping implements the Neumann condition at second order too — it is not degrading the
boundary rows. Doubling κ doubles the error to 2.000000, as it must for a linear operator.

Probe: `analysis/probes/probeMMSCahnHilliard.m`; data `analysis/probes/data/mms_cahnhilliard.{mat,csv}`.

**Gate 3 (κ contrast in a verified-active regime) — passed.** 3 d, N = 30. The probe reports
its own regime before drawing any conclusion: `phib_sup(κ=1e-7) = 0.1419`,
`mat_mass = 0.0363 kg/m²` — the supernatant genuinely carries biofilm, so the contrast is
meaningful.

| quantity | κ = 1e-7 | κ = 1e-2 | rel. diff |
|---|---|---|---|
| `phib_sup` | 0.141885 | 0.048612 | 0.657 |
| `mat_mass` | 0.0363 | 0.0668 | 0.840 |
| `mat_depth` | 0.0328 | 0 | 1.000 |

Direction is physically sensible: stronger surface tension smears the mat, dropping the peak
volume fraction (so `mat_depth`, a `φ_b > 0.05` band, vanishes) while retaining more total
biomass.

**Gate 4 (dt cost) — measured, and the prediction was wrong.** With `AdaptiveMaxDt` at its
usual 3e-6 the step counts were identical (1065422 vs 1065432) — but that only measures the
cap: 3 d / 3e-6 = 1.0e6 steps, so dt was cap-bound throughout and the CFL never bound.
Re-run with the cap lifted to 1e-3 so the CFL binds:

| N | κ = 1e-7 | κ = 1e-2 | ratio |
|---|---|---|---|
| 30 | 113141 | 113359 | 1.002 |
| 100 | 133367 | 133919 | 1.004 |

**0.2–0.4 %**, even at κ = 1e-2 (10⁵× physical). The plan predicted ~2×, estimating the new
`v_b` term as `ζ₀κ/Δz³` on the assumption that `u` varies over one cell. It does not: the
interface is unresolved, so the numerical front is smooth over many cells and the third
derivative of `u` is far smaller than `1/Δz³`. **The dz⁻⁴ penalty is latent** and will appear
exactly when a mesh resolves the interface — i.e. in the graded-mesh work parked in
`PLANS.md`, whose cost estimate should be revised upward accordingly.

**Julia mirror — done, and agrees.** `SSF.jl/src/SandFilter.jl`: same `+ n0` row shift, same
`gridsize`/`gridzero` fixes plus `cellwidths`/`centerspacings`. Cross-checked against MATLAB
on the same configuration: `n0 = 21`, block (2,1) `nnz = 61`, block (2,2) `== I`, block (1,1)
free of κ — identical in both.

**Goldens — all four pass unchanged, and that is itself a finding.** They were expected to
break. They do not, because they run a near-uniform state over `SimulationTime = 1e-3`, where
the curvature `L·u ≈ 0` and the κ term is inert *whether or not it is wired up*. **The golden
suites are structurally incapable of detecting this class of defect** — the same
inactive-regime trap that invalidated the first test above. A regression test that does catch
it has been added (`test/runtests.jl`, "Cahn-Hilliard: kappa actually reaches mu"): it forces
a non-uniform `u`, asserts `μ ≠ Ψ'(u)`, and asserts the departure scales linearly in κ.

**The old unit test encoded the defect.** `test/runtests.jl:165` asserted
`all(mats.diffusion.I .<= n0)` under the comment "Convection and diffusion live in the first
block". That is the wrong placement stated as an invariant, which is presumably how the
misplacement survived. Updated.

## Goldens re-exported (2026-08-23)

All four reference suites in `SSF.jl/test/golden/` were regenerated after both MATLAB
trees were fixed — `reference`/`reference_adaptive` from MPC-SSF, `reference_pathogen`/
`reference_pathogen_light` from the legacy tree via `run_export_pathogen.m`. All exports
returned `flag = OK`. Numbers, tables and the full method are in
`SSF.jl/test/golden/README.md` under "Re-export 2026-08-23"; the headline:

- `reference_adaptive` came back **bit-identical** (its spec runs 1e-5 d).
- The pathogen references moved in `vel_biofilm` (maxabs 1.07e-05, rel 1.5e-06) — the
  expected signature, since `v_b` is Solver A's only output and κ enters it through μ.
- Julia–MATLAB agreement on the fixed-step suite went from maxabs **8.946e-09** to
  **8.882e-16**. Julia already carried the fix, so that 8.9e-09 *was* the κ term; with
  both sides fixed it closes to machine precision. Adaptive and pathogen residuals are
  unchanged, being dominated by pre-existing roundoff/cross-codebase differences.
- `Pkg.test`: 285/285 pass, all four suites `RESULT: MATCH`.

This also confirms the earlier reading that **the goldens cannot discriminate κ** — they
pass both before and after because their short, near-uniform runs leave the term small
either way. Discrimination is `test/cohesion_discrimination.jl` and the
`Cahn-Hilliard: kappa actually reaches mu` testset.

## Consequences

- **All four golden suites must be re-exported.** They pass today only because κ is inert.
- **Timestep.** κ re-enters the CFL through `vbmax`: `v_b` gains `~ζ₀κ/Δz³`, so
  `w_v/Δz ~ ζ₀κ/Δz⁴`. Measured, not assumed — see the run log.
- **Physics.** Surface tension that never acted now acts, in the supernatant only
  (`v_b ≡ 0` for z > 0). Expect movement in `phib_sup`, `mat_mass`, `mat_depth`; bed
  quantities should be far less affected.
- **Everything reported to date was produced with κ inert.** The cover-sweep schmutzdecke
  results (`R_mat` 0.25 at Stage A, 0.001 at Stage B) sit exactly in the supernatant, where
  the term acts, and need re-running. The earlier concern that they were mesh-artefacts was
  the right worry for the wrong reason: Appendix A.2's Δz < √κ criterion was not being
  *violated*, it was **moot**, because no interfacial layer was being computed.

## Alternatives rejected

- **Reinstate the Dirichlet boundary term as well.** Rejected: Appendix A.1 records the
  switch to Neumann as deliberate and motivated by timestep cost. Reinstating it would
  silently revert a published modelling decision.
- **Rename `GridSize` throughout.** Rejected: ~60 call sites, all correct on a uniform mesh.
  Making it error on non-uniform meshes gives the same protection with no churn.
- **Replace the scheme with Bailo2023 now.** Deferred to `PLANS.md`: the right first step is
  to make the existing scheme do what it was documented to do, so the change in behaviour
  attributable to κ can be measured against a known baseline.
