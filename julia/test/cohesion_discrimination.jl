# Tests that discriminate the cohesion potential.
#
# WHY THIS FILE EXISTS. The Julia port ran for months with a generic double well,
# 0.25 u^2 (1-u)^2, instead of the published dpsi/du = u^2 (u - 3 zeta_1 / 2) —
# a discrepancy reaching 98x at low phi and inverting near the self-limiting peak.
# All four golden-master suites passed unchanged across that change, because two
# start from clean initial conditions where phi_b = 0 (both forms vanish there)
# and the pathogen goldens run only 1e-5 days. Nothing in the suite could see it.
#
# The bug was never in the formula — it was in the wiring. The published form
# existed in `CahnHilliardModel`, and no preset carried it into the solver. So
# these tests check the wiring, not the algebra:
#
#   1. Every preset's potential IS the published form, evaluated at that preset's
#      own zeta_1, and is NOT the double well. This is the direct guard on the
#      defect that actually occurred.
#   2. zeta_1 regenerates the handle when a model is reconstructed correctly, and
#      the documented stale-closure trap is pinned so it stays visible: because
#      `potential_gradient` closes over `zeta_1`, a caller that reconstructs while
#      passing the old handle changes the stored number and not the behaviour.
#
# Together with the `CahnHilliardModel` unit tests in runtests.jl, which pin the
# closed form itself, a regression in the potential or its wiring fails loudly.
#
# NOT COVERED HERE: a solver-level assertion that the potential changes simulated
# output. See "Solver-level discrimination" at the foot of this file — the regime
# it requires is established and measured, but it needs a mature-biofilm state
# that is too slow to build inside the unit suite.

using MPCSSF
using Test

# The published potential, written out independently of the implementation so
# these tests cannot pass by construction if the source drifts.
_published_dpsi(u, zeta_1) = u^2 * (u - 1.5 * zeta_1)

# The wrong form the port used, kept as the adversary.
_double_well_dpsi(u) = 0.25 * u^2 * (1 - u)^2

@testset "cohesion potential is discriminated" begin

    @testset "presets carry the published potential" begin
        for (name, m) in (("simpleModel", simpleModel()),
                          ("modelLund", modelLund()),
                          ("modelPathogen", modelPathogen()))
            ch = m.cohesion_submodel
            @test ch !== nothing
            z1 = ch.zeta_1
            for u in (0.01, 0.05, 0.30, 0.63, 0.90)
                @test ch.potential_gradient(u) ≈ _published_dpsi(u, z1) rtol=1e-12
            end
            # phi = 0.01 is where the two forms differ most (98x), so it is the
            # sharpest single check that this is not the double well.
            @test !isapprox(ch.potential_gradient(0.01), _double_well_dpsi(0.01);
                            rtol=1e-6)
        end
    end

    @testset "zeta_1 regenerates the handle" begin
        # Reconstructed WITHOUT passing potential_gradient, which is the documented
        # correct usage: the handle regenerates against the new zeta_1.
        lo = CahnHilliardModel(kappa=1e-2, zeta_0=1.0, zeta_1=0.001)
        hi = CahnHilliardModel(kappa=1e-2, zeta_0=1.0, zeta_1=0.05)
        @test lo.potential_gradient(0.05) ≈ _published_dpsi(0.05, 0.001) rtol=1e-12
        @test hi.potential_gradient(0.05) ≈ _published_dpsi(0.05, 0.05) rtol=1e-12
        @test lo.potential_gradient(0.05) != hi.potential_gradient(0.05)

        # The trap, pinned deliberately: reconstructing while carrying the OLD
        # handle over changes the stored zeta_1 and nothing else. Asserted so the
        # failure mode stays visible rather than becoming a surprise later.
        stale = CahnHilliardModel(kappa=1e-2, zeta_0=1.0, zeta_1=0.05,
                                  potential_gradient=lo.potential_gradient)
        @test stale.zeta_1 == 0.05
        @test stale.potential_gradient(0.05) == lo.potential_gradient(0.05)
        @test stale.potential_gradient(0.05) != hi.potential_gradient(0.05)
    end
end

# ── Solver-level discrimination: measured, not yet automated ─────────────────
#
# Whether the potential reaches the solver at all can only be shown by running
# one. Measured 2026-07-29, comparing the published form against the double well
# on otherwise identical runs, as max|dphi_b| / max(phi_b) in the final frame:
#
#   initial condition            phi_b      relative difference
#   ---------------------------  ---------  -------------------
#   clean start (SimpleModel)    0          0          (both forms vanish)
#   seeded bump, amp 5e-3        4.9e-3     7e-17      (machine epsilon)
#   same, zeta_0 raised to 1e6   4.9e-3     1e-10
#   MATURE biofilm, 1 day        1.6e-1     1.0e-1 to 1.3e-1
#
# The cohesion flux scales with the potential and with mobility u(1-u), so at
# phi_b ~ 5e-3 it is negligible no matter how large zeta_0 is: raising zeta_0 by
# 1e6 moves the difference from 1e-17 only to 1e-10. A mature biofilm at
# phi_b ~ 0.16 discriminates by 10-13%, nine orders of magnitude better.
#
# This is the quantitative reason all four goldens were blind, and it says what a
# solver-level test needs: a mature biofilm, not a seeded one. Building that
# state means a growth pre-run of order a simulated day, which is too slow for
# the unit suite. The right home is a golden with a committed mature initial
# state, anchored to MATLAB — the legacy implementation was always correct, so a
# reference can be exported from it. That is the outstanding piece of work.
