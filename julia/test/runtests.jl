using MPCSSF
using Test
using SparseArrays   # for `sparse(::Triplets)` in the Cahn-Hilliard tests

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

    @testset "CahnHilliardModel" begin
        ch = CahnHilliardModel(kappa=1e-3, zeta_0=2.0)
        @test ch.kappa == 1e-3
        @test ch.zeta_1 == 0.0                  # default
        # Default mobility u(1-u) and potential gradient, scalar + vector.
        @test ch.mobility(0.5) ≈ 0.25
        @test ch.mobility([0.0, 0.5, 1.0]) ≈ [0.0, 0.25, 0.0]
        @test ch.potential_gradient(0.5) ≈ 0.25 * (0.25 * 0.25)
    end

    @testset "SandFilter grid" begin
        f = SandFilter(height=1.0, depth=1.0)
        @test f.grid === nothing
        f = addgridpoints(f, 50)
        @test f.grid isa Grid
        b = f.grid.boundaries
        @test issorted(b)                                   # ascending faces
        @test length(f.grid.centers) == length(b) - 1
        @test f.grid.centers ≈ 0.5 .* (b[1:end-1] .+ b[2:end])
        @test b[1] ≈ -f.height                              # spans up from −height
        # n0 marks the cell straddling z = 0.
        n0 = gridzero(f)
        @test abs(f.grid.centers[n0]) < gridsize(f) / 2
    end

    @testset "computeporosity" begin
        f = SandFilter(sand_porosity=0.4, sand_roughness=5e-3)
        @test computeporosity(f, 0.0) ≈ 0.4          # sand surface
        @test computeporosity(f, 1.0) ≈ 0.4          # deep in sand (clamped)
        @test computeporosity(f, -1.0) ≈ 1.0         # supernatant water (clamped)
        @test computeporosity(f, -0.002) ≈ 0.64      # within roughness ramp
        @test computeporosity(f, [0.0, -1.0]) ≈ [0.4, 1.0]   # vectorized
    end

    @testset "Cahn-Hilliard matrices" begin
        f = addgridpoints(SandFilter(), 40)
        model = Model(cohesion_submodel=CahnHilliardModel(kappa=1e-3, zeta_0=2.0))
        mats = get_cahn_hilliard_matrices(f, model, false)
        n0 = gridzero(f)
        N = 2n0
        for t in (mats.convection, mats.diffusion, mats.diffusion_mobility)
            @test (t.m, t.n) == (N, N)
            S = sparse(t)
            @test size(S) == (N, N)
        end
        # Convection and diffusion live in the first block (rows & cols ≤ n0).
        @test all(mats.convection.I .<= n0) && all(mats.convection.J .<= n0)
        @test all(mats.diffusion.I .<= n0) && all(mats.diffusion.J .<= n0)
        # Mobility-weighted diffusion couples block 1 (rows) to block 2 (cols > n0).
        @test all(mats.diffusion_mobility.I .<= n0)
        @test all(mats.diffusion_mobility.J .> n0)
        # Upwind vs. centered convection differ.
        up = get_cahn_hilliard_matrices(f, model, true)
        @test up.convection.V != mats.convection.V
    end

    @testset "Model" begin
        het = Particle(name="HET", density=1.0)
        pom = Particle(name="POM", density=2.0)
        dom = Liquid(name="DOM", density=998.0)
        rx1 = Reaction(name="growth", nominal_rate=2.0,
                       temperature_correction_factor=1.07,
                       stoichiometric_coefficients=Dict("HET"=>1.0, "DOM"=>-0.5),
                       half_saturation_constants=Dict("DOM"=>0.3))
        rx2 = Reaction(name="decay", nominal_rate=0.5,
                       stoichiometric_coefficients=Dict("HET"=>-1.0, "POM"=>1.0))

        # Positional constructor + defaults.
        m = Model([het, pom, dom], [rx1, rx2])
        @test m.water_density == 998.0
        @test m.biofilm_porosity == 0.99
        @test m.osmosis_rate == 1e-5
        @test m.cohesion_submodel === nothing

        # Component views preserve declaration order.
        @test particles(m) == [het, pom]
        @test liquids(m) == [dom]
        @test eltype(particles(m)) == Particle

        # Reaction rates at the nominal temperature reduce to nominal rates.
        @test compute_reaction_rates(m, 20) ≈ [2.0, 0.5]

        # Derived matrices have the right shapes and agree with the lookups.
        @test size(stoichiometric_matrix_particles(m)) == (2, 2)
        @test size(stoichiometric_matrix_liquids(m)) == (1, 2)
        @test size(half_saturation_constants(m)) == (3, 2)
        @test stoichiometric_matrix_particles(m) ==
              lookup_stoichiometric_coefficients(m.reactions, particles(m))
        # HET stoichiometry across the two reactions: +1 (growth), −1 (decay).
        @test stoichiometric_matrix_particles(m)[1, :] == [1.0, -1.0]
        # DOM (the only liquid) is consumed in reaction 1.
        @test stoichiometric_matrix_liquids(m)[1, :] == [-0.5, 0.0]

        # Empty model is well-defined.
        @test isempty(particles(Model()))
        @test compute_reaction_rates(Model(), 20) == Float64[]
    end

    @testset "State" begin
        f = addgridpoints(SandFilter(), 30)
        het = Particle(name="HET", density=1100.0)
        pom = Particle(name="POM", density=1100.0)
        dom = Liquid(name="DOM", density=998.0)
        m = Model([het, pom, dom])
        s = State(f, m)

        N = length(f.grid.centers)
        kP, kL = 2, 1
        @test s.time == 0.0
        @test size(s.global_concentration.matrix) == (N, kP)
        @test size(s.global_concentration.enclosed_liquids) == (N, kL)
        @test length(s.velocity.biofilm) == N - 1
        @test length(s.velocity.flowing) == N + 1
        @test length(s.enclosed_water_volume) == N

        @test size(global_concentration_biofilm(s)) == (N, 2kP + kL)
        @test size(global_concentration_flowing(s)) == (N, kP + kL)
        @test all(==(0), global_concentration_biofilm(s))      # clean filter

        # Volume fraction = concentration / density (mutate a cell to check).
        s.global_concentration.matrix[1, 1] = 1100.0
        @test volume_fractions(s).matrix[1, 1] ≈ 1.0
        @test size(volume_fractions(s).enclosed_liquids) == (N, kL)

        # A filter without a grid cannot form a State.
        @test_throws Exception State(SandFilter(), m)
    end

    @testset "simulate (to port)" begin
        f = addgridpoints(SandFilter(), 20)
        m = Model([Particle(name="HET", density=1.0)])
        s = State(f, m)
        @test_throws Exception simulate(s)   # solver not yet ported
    end
end
