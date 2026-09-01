# Grid refinement drivers + parallel Dirichlet-BC study

## Goal

Two coupled questions, one driver family.

1. **Refinement.** Does anything converge under dz -> 0 near z = 0? The 0-cell sits
   exactly on the point where v_b is undefined (model.tex:148-152), so phi_b(0) is
   not obviously a convergent quantity. Measure which outputs converge and which do
   not.
2. **Boundary condition.** The current z = 0 condition is homogeneous Neumann
   d(phi_b)/dz = 0 (index clamping in getCahnHilliardMatrices.m, manuscript
   Appendix A.1). The superseded Diehl2025 variant imposes phi_b continuity by
   taking the ghost value from the first BED cell, which follows its own ODE and is
   therefore GIVEN. Run both and compare.

Done = a convergence table over N for both BCs, and a defensible statement of which
quantities are grid-independent.

## Key facts established before writing any code

- Grid: dz = Depth/(N + 1/2) puts a cell CENTRE at z = 0 and no face there
  (SandFilter.m:119-127). n0 = index of that cell.
- Current BC: the `min(max(1,...),n0)` clamping makes the k=0 and k=n0 faces
  contribute +1-1-1+1 = 0, i.e. homogeneous Neumann at BOTH ends.
- Dirichlet variant (derived from the code's own convention, and matching
  numerical-method.md:250-262):
    row n0 of the mu block gains -kappa/dz^2 on the diagonal, and
    rhs(mu row n0) gains -kappa*u_bed/dz^2,   u_bed = phi_b of the first bed cell.
  Setting u_bed := u(n0) must reduce EXACTLY to the Neumann case. That identity is
  the self-test, so the sign does not rest on my algebra.
- NO spike exists today: phi_b(0-cell)/phi_b(neighbour) = 1.011 at N=100, 20 d,
  because cohesion is in its spreading branch (Psi'' = +0.298 at phi = 0.32,
  D_eff ~ 6.5 m^2/d, an 11 m diffusion length in a 1 m supernatant). The
  surface-source concern is DORMANT, not absent -- it can wake when zeta_1 rises or
  when the BC changes.

## Steps

1. **`analysis/probes/probeRefinement.m`** -- one (N, BC) pair per call; saves the
   profile plus derived metrics. No fixed length anywhere in the metrics:
   - Sigma = sum(eps.*phi_b)*dz over the supernatant  [convergent by construction]
   - phiPeak, and equivalent thickness L = Sigma/phiPeak
   - front position where phi_b crosses zeta_1 (self-consistent: that is where the
     potential changes character)
   - sharpness = width over which phi_b falls 0.9*peak -> 0.1*peak
   - zeroCellRatio = phi_b(n0)/phi_b(n0-1)   <-- the delta detector
   - bed penetration depth (where bed phi_b falls to 10% of its top value)
   *Verify:* checkcode clean; metrics reproduce the known N=100 values
   (Sigma = 0.0839, ratio = 1.011) when run at N=100, 20 d.

2. **`CohesionBC` parameter** in simulate.m, "neumann" (default, current behaviour)
   or "dirichlet". Implemented as a one-off modification of CH0 plus a per-step rhs
   term; nothing else in the hot loop changes.
   *Verify:* the self-test above -- `BCSelfTest=true` forces u_bed := u(n0) and the
   run must match "neumann" BIT-FOR-BIT. If it does not, the sign is wrong.

3. **`slurm/refinement.sbatch`** -- array over (N, BC), N in {100, 200, 500, 1000}.
   *Verify:* all arms Flag == "OK"; step counts scale as expected once CFL-bound.

4. **Read the table.** Converged: Sigma, L, front position, deposition flux.
   Suspect: phi_b(0-cell), and any reaction rate evaluated there.
   *Verify:* zeroCellRatio vs N. Flat => no delta forming, refinement is meaningful.
   Growing like 1/dz => the surface source is real and phi_b(0) has no limit.

## Files touched

- `analysis/probes/probeRefinement.m` (new)
- `src/@State/simulate.m` (CohesionBC, BCSelfTest; diagnostic-free hot path)
- `slurm/refinement.sbatch` (new)

## Deliberately NOT done

- No grid change to put a face at z = 0. The Diehl variant does not need one: the
  ghost is the bed cell, which already exists.
- No change to the default BC. Neumann stays default; Dirichlet is opt-in, exactly
  as getCahnHilliardMatrices.m:57-61 warns ("must NOT be reinstated here").
