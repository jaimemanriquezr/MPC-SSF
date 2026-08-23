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
