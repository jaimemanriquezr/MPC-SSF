using MPCSSF
using Test

# Phase 1 verification scaffold. As each unit is ported, the to-port testsets
# are filled in. Ecological layer is implemented; the solver remains a stub.
# Where possible, replace these with golden-master comparisons against MATLAB
# reference outputs for SimpleModel / LundMPCModel (outline §5).

@testset "MPCSSF" begin

    @testset "components" begin
        p = Particle(name="HET", density=1.0, attachment_sand=2.0)
        @test p isa Component
        @test p.name == "HET"
        @test p.density == 1.0
        @test p.attachment_sand == 2.0
        @test p.transport_rate == 0.0        # default

        l = Liquid(name="DOM", density=1000.0, dispersivity=1e-3)
        @test l isa Component
        @test l.dispersivity == 1e-3
    end

    @testset "compute_rate" begin
        rx = Reaction(nominal_rate=2.0, temperature_correction_factor=1.07)
        # At the nominal temperature (20 °C) the correction is θ^0 = 1.
        @test compute_rate(rx, 20) ≈ 2.0
        # Warmer ⇒ faster (θ > 1); cooler ⇒ slower.
        @test compute_rate(rx, 30) > 2.0
        @test compute_rate(rx, 10) < 2.0
        # Kelvin scale with matching nominal reproduces the nominal rate.
        @test compute_rate(rx, 293; scale=:kelvin, nominal_temperature=293) ≈ 2.0
        @test_throws ArgumentError compute_rate(rx, 20; scale=:fahrenheit)
        # Broadcasts over a vector of reactions.
        rxs = [rx, Reaction(nominal_rate=5.0)]
        @test compute_rate.(rxs, 20) ≈ [2.0, 5.0]
    end

    @testset "lookup matrices" begin
        het = Particle(name="HET", density=1.0)
        pom = Particle(name="POM", density=2.0)
        dom = Liquid(name="DOM", density=1000.0)
        comps = Component[het, pom, dom]

        # 2 HET -> POM, consuming DOM; rate ~ HET (order 1), Monod in DOM.
        rx = Reaction(name="growth",
                      stoichiometric_coefficients=Dict("HET"=>-1.0, "POM"=>0.5),
                      order=Dict("HET"=>1.0),
                      half_saturation_constants=Dict("DOM"=>0.3))
        rxs = [rx]

        sigma = lookup_stoichiometric_coefficients(rxs, comps)
        @test size(sigma) == (3, 1)
        @test sigma[:, 1] == [-1.0, 0.5, 0.0]      # absent ⇒ 0

        p = lookup_order(rxs, comps)
        @test p[:, 1] == [1.0, 0.0, 0.0]

        K = lookup_half_saturation_constants(rxs, comps)
        @test K[3, 1] == 0.3
        @test isnan(K[1, 1]) && isnan(K[2, 1])     # absent ⇒ NaN

        # Component-keyed dicts are accepted (keys normalized to names).
        rx2 = Reaction(stoichiometric_coefficients=Dict(het=>-1.0, pom=>0.5))
        @test lookup_stoichiometric_coefficients([rx2], comps)[:, 1] == [-1.0, 0.5, 0.0]
    end

    @testset "lookup_quotients" begin
        het = Particle(name="HET", density=1.0)
        pom = Particle(name="POM", density=2.0)
        dom = Liquid(name="DOM", density=1000.0)
        comps = Component[het, pom, dom]

        # An inhibition/quotient term DOM/HET in a single reaction.
        rx = Reaction(half_saturation_constants=Dict("DOM/HET"=>0.7))
        q = lookup_quotients([rx], comps)
        @test size(q.K) == (1, 1)
        @test q.K[1, 1] == 0.7
        @test q.num_idx == [3]   # DOM is component 3
        @test q.den_idx == [1]   # HET is component 1

        # No quotient keys ⇒ empty result.
        q0 = lookup_quotients([Reaction(half_saturation_constants=Dict("DOM"=>0.3))], comps)
        @test size(q0.K) == (0, 1)
        @test isempty(q0.num_idx)
    end

    @testset "solver (to port)" begin
        f = SandFilter()
        m = Model()
        @test_throws Exception State(f, m)   # constructor not yet ported
    end
end
