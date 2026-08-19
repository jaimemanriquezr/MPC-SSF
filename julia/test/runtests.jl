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
        # Published potential dpsi/du = u^2 (u - 3 zeta_1 / 2); degenerates to u^3
        # when zeta_1 = 0 (see .claude/decisions/2026-07-28-published-cohesion-potential.md).
        @test ch.potential_gradient(0.5) ≈ 0.125
        @test ch.potential_gradient([0.0, 0.5]) ≈ [0.0, 0.125]
        # With zeta_1 set it must reproduce the handles baked into the legacy .mat
        # models, and zeta_1 must actually move the potential.
        ch2 = CahnHilliardModel(kappa=1e-7, zeta_0=1e3, zeta_1=0.005)
        @test ch2.potential_gradient(0.05) ≈ 0.05^2 * (0.05 - 3 * (1 / 200) / 2)
        @test ch2.potential_gradient(0.63) ≈ 0.63^2 * (0.63 - 0.0075)
        @test ch2.potential_gradient(0.05) != ch.potential_gradient(0.05)
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

    @testset "simulate end-to-end (smoke)" begin
        # A zero-rate model on a clean filter: exercises the whole solver path
        # (matrix assembly, SOLVER A linear solve, SOLVER B transport, frame
        # capture) while staying numerically trivial (no sources ⇒ stays clean).
        f = addgridpoints(SandFilter(), 20)
        het = Particle(name="HET", density=1000.0)
        dom = Liquid(name="DOM", density=998.0)
        rx = Reaction(name="noop", nominal_rate=0.0, optimal_light_factor=1.0,
                      order=Dict("HET"=>1.0),
                      stoichiometric_coefficients=Dict("HET"=>1.0, "DOM"=>-1.0))
        ch = CahnHilliardModel(kappa=1e-3, zeta_0=0.5)
        m = Model([het, dom], [rx]; cohesion_submodel=ch)
        s = State(f, m)

        res = simulate(s; simulation_time=1e-4, time_step=1e-5, n_frames=5, quiet=true)
        N = length(f.grid.centers)
        @test res.flag == "OK"
        @test size(res.frames[:concentration_biofilm]) == (N, 5, 3)   # 2kP+kL
        @test size(res.frames[:concentration_flowing]) == (N, 5, 2)   # kP+kL
        @test length(res.frames[:time]) == 5
        @test all(==(0), res.frames[:concentration_biofilm])          # clean stays clean
        @test haskey(res.simulation_data, :time_final)
        # Biofilm velocity in the supernatant is the advective velocity (no CH forcing).
        @test any(!=(0), res.frames[:velocity_biofilm])
    end

    @testset "Results accessors" begin
        f = addgridpoints(SandFilter(), 20)
        het = Particle(name="HET", density=1000.0)
        dom = Liquid(name="DOM", density=998.0)
        rx = Reaction(name="noop", nominal_rate=0.0, optimal_light_factor=1.0,
                      order=Dict("HET"=>1.0))
        ch = CahnHilliardModel(kappa=1e-3, zeta_0=0.5)
        m = Model([het, dom], [rx]; cohesion_submodel=ch)
        res = simulate(State(f, m); simulation_time=1e-4, time_step=1e-5,
                       n_frames=4, quiet=true)

        N = length(f.grid.centers)
        @test length(times(res)) == 4
        @test depths(res) == f.grid.centers
        @test size(concentration(res, "HET", :matrix)) == (N, 4)
        @test size(concentration(res, "HET", :flowing)) == (N, 4)
        @test all(==(0), concentration(res, "DOM", :matrix))      # liquids: no matrix
        @test size(concentration(res, "Water", :enclosed)) == (N, 4)
        @test_throws Exception concentration(res, "HET", :nope)
        @test_throws Exception concentration(res, "Ghost", :flowing)

        vf = get_volume_fractions(res)
        @test size(vf.biofilm) == (N, 4)
        @test all(==(0), vf.biofilm)                              # clean filter
        @test vf.biofilm == vf.matrix .+ vf.enclosed

        # Set HET matrix concentration to its density ⇒ φ_M = 1 in that cell/frame.
        res.frames[:concentration_biofilm][1, 1, 1] = 1000.0
        @test get_volume_fractions(res).matrix[1, 1] ≈ 1.0
    end

    @testset "plotting" begin
        # Reaction-rate frame accessors (feed plot_reaction_rates). Two named
        # reactions so column order and labels can be checked.
        f = addgridpoints(SandFilter(), 20)
        het = Particle(name="HET", density=1000.0)
        dom = Liquid(name="DOM", density=998.0)
        rx1 = Reaction(name="growth", nominal_rate=0.0, optimal_light_factor=1.0,
                       order=Dict("HET"=>1.0))
        rx2 = Reaction(name="decay", nominal_rate=0.0, order=Dict("HET"=>1.0),
                       stoichiometric_coefficients=Dict("HET"=>-1.0))
        ch = CahnHilliardModel(kappa=1e-3, zeta_0=0.5)
        m = Model([het, dom], [rx1, rx2]; cohesion_submodel=ch)
        res = simulate(State(f, m); simulation_time=1e-4, time_step=1e-5,
                       n_frames=4, quiet=true)

        N = length(f.grid.centers)
        rr = reaction_rates(res)
        @test size(rr) == (N, 4, 2)                  # depth × frame × reaction
        @test rr === res.frames[:reaction_rates]
        @test all(==(0), rr)                          # zero nominal rate ⇒ no reaction
        @test reaction_names(res) == ["growth", "decay"]

        # Without a Makie backend loaded, every plot generic hits the friendly
        # fallback: an ArgumentError-free `error` naming the function and Makie.
        for (fn, call) in (
                ("plot_concentration",         () -> plot_concentration(res, "HET", :flowing)),
                ("plot_concentration_heatmap", () -> plot_concentration_heatmap(res, "HET", :flowing)),
                ("plot_volume_fractions",      () -> plot_volume_fractions(res)),
                ("plot_velocity",              () -> plot_velocity(res)),
                ("plot_cfl",                   () -> plot_cfl(res)),
                ("plot_reaction_rates",        () -> plot_reaction_rates(res)))
            err = try; call(); nothing; catch e; e; end
            @test err isa ErrorException
            @test occursin(fn, err.msg) && occursin("Makie", err.msg)
        end
    end

    @testset "modelLund preset" begin
        m = modelLund()
        # 9 components (4 particulate, 5 dissolved), 5 reactions.
        @test length(m.components) == 9
        @test length(particles(m)) == 4
        @test length(liquids(m)) == 5
        @test [c.name for c in particles(m)] == ["HET", "PHO", "POM", "PAT"]
        @test [c.name for c in liquids(m)] == ["O2", "IC", "NH4", "HPO4", "DOM"]
        @test length(m.reactions) == 5

        # Global params + cohesion submodel (verbatim from modelLund.m).
        @test m.biofilm_porosity == 0.99
        @test m.osmosis_rate == 1.00e-7
        @test m.cohesion_submodel.kappa == 1.00e-6
        @test m.cohesion_submodel.zeta_0 == 1.00e6
        @test m.cohesion_submodel.zeta_1 == 1/100
        # Detachment @(v) sqrt(v/7.2): scalar and vector.
        @test m.detachment(7.2) ≈ 1.0
        @test m.detachment([7.2, 28.8]) ≈ [1.0, 2.0]

        # Particle physical params (shared across the four particulates).
        het = particles(m)[1]
        @test het.density ≈ 1.117e3
        @test het.attenuation ≈ 0.094
        @test het.transport_rate ≈ 5.47
        @test het.attachment_sand ≈ 5.47e2
        # POM does not attach.
        @test particles(m)[3].attachment_sand == 0.0

        # Derived matrices have the right shape and a couple of known entries.
        σ = stoichiometric_coefficients(m)          # 9 × 5
        @test size(σ) == (9, 5)
        @test σ[1, 1] == 1.0                         # HET in heterotroph growth
        @test σ[5, 1] ≈ -1.2317                      # O2 in heterotroph growth
        @test size(stoichiometric_matrix_particles(m)) == (4, 5)
        @test size(stoichiometric_matrix_liquids(m)) == (5, 5)
        @test size(reaction_orders(m)) == (9, 5)

        # Hydrolysis carries a POM/HET quotient (5th reaction).
        q = quotients(m)
        @test length(q.num_idx) == 1
        @test q.num_idx[1] == 3                       # POM
        @test q.den_idx[1] == 1                       # HET
        @test q.K[1, 5] ≈ 2.00e-5

        # Phototroph growth is the only light-dependent reaction.
        @test [r.is_light_dependent for r in m.reactions] == [false, true, false, false, false]
    end

    @testset "phototroph respiration (metabolic split)" begin
        # Off by default: modelLund() is the original 5-reaction model.
        m0 = modelLund()
        @test length(m0.reactions) == 5
        @test m0.reactions[2].minimum_light_factor == 0.01

        # On: sixth reaction, reversed stoichiometry, dark floor retired.
        m = modelLund(phototroph_respiration=0.276)
        @test length(m.reactions) == 6
        resp = m.reactions[6]
        @test resp.name == "Phototroph respiration"
        @test resp.nominal_rate == 0.276
        @test !resp.is_light_dependent
        @test resp.light_inhibition ≈ 8e-5/1.814e-2   # Wolf2007 K_inh, normalized
        @test_throws ArgumentError Reaction(name="bad", is_light_dependent=true,
                                            light_inhibition=1e-3)
        gro = m.reactions[2]
        @test gro.minimum_light_factor == 0.0
        for k in ("O2", "IC", "NH4", "HPO4")
            @test resp.stoichiometric_coefficients[k] ==
                  -gro.stoichiometric_coefficients[k]
        end
        # Dark column: with respiration ON, oxygen is CONSUMED where algae sit;
        # with it OFF (floor 0) the O2 field is untouched by phototrophs. Seed
        # PHO + O2 via the influent, kill the light, run briefly, compare total
        # flowing+enclosed O2.
        infl = [0.0, 1.0e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.0e-5, 0.0, 0.0]
        function darkrun(m)
            f = addgridpoints(SandFilter(), 15)
            f.light_irradiation = t -> 0.0
            simulate(State(f, m); inflow_concentrations=infl, simulation_time=2e-3,
                     time_step=:adaptive, adaptive_initial_dt=1e-8, adaptive_max_dt=3e-6,
                     n_frames=3, implicit_osmosis=true, quiet=true)
        end
        r_on  = darkrun(m)
        r_off = darkrun(modelLund(phototroph_respiration=0.0))
        o2 = r -> sum(concentration(r, "O2", :flowing)[:, end]) +
                  sum(concentration(r, "O2", :enclosed)[:, end])
        @test r_on.flag == "OK" && r_off.flag == "OK"
        @test o2(r_on) < o2(r_off)          # respiration consumes O2 in the dark
        # PHO mass: respiration consumes biomass while the floor model's dark
        # growth adds it — both effects separate the runs in the same direction.
        # (IC is not asserted: over this short horizon the advected influent IC
        # swamps the reaction signal.)
        pho = r -> sum(concentration(r, "PHO", :matrix)[:, end]) +
                   sum(concentration(r, "PHO", :enclosed)[:, end]) +
                   sum(concentration(r, "PHO", :flowing)[:, end])
        @test pho(r_on) < pho(r_off)        # respiration removes biomass

        # Dark switch (Wolf2007 r6): under bright constant light the inhibition
        # factor K/(K+I) collapses (K = 4.4e-3 vs I ~ 1), so a respiration-only
        # model barely touches O2; in darkness it consumes it.
        pho_c = Particle(name="PHO", density=1.117e3, dispersivity=1.2e-2, transport_rate=5.47)
        o2_c  = Liquid(name="O2", density=998.0, dispersivity=1.2e-2, transport_rate=600.0)
        resp_only = Reaction(name="resp", nominal_rate=0.55, order=Dict("PHO"=>1.0),
                             light_inhibition=8e-5/1.814e-2,
                             stoichiometric_coefficients=Dict("PHO"=>-1.0, "O2"=>-0.9301))
        mR = Model([pho_c, o2_c], [resp_only];
                   cohesion_submodel=CahnHilliardModel(kappa=1e-6, zeta_0=1e2, zeta_1=1e-2))
        function lightrun(lum)
            f = addgridpoints(SandFilter(), 15)
            f.light_irradiation = t -> lum
            # Fixed step: the adaptive-CFL bound assumes the Lund reaction
            # structure and cannot run this 2-component model.
            simulate(State(f, mR); inflow_concentrations=[1e-2, 9.1e-3],
                     simulation_time=2e-3, time_step=3e-6, n_frames=3,
                     implicit_osmosis=true, quiet=true)
        end
        o2tot = r -> sum(concentration(r, "O2", :flowing)[:, end]) +
                     sum(concentration(r, "O2", :enclosed)[:, end])
        rDark = lightrun(0.0); rBright = lightrun(1.0)
        @test rDark.flag == "OK" && rBright.flag == "OK"
        drop_dark = 1 - o2tot(rDark)/o2tot(rBright)
        @test o2tot(rDark) < o2tot(rBright)          # dark consumes more
        @test drop_dark > 0                          # sanity on the sign
    end

    @testset "modelPathogen preset" begin
        m = modelPathogen()
        # 9 components (4 particulate, 5 dissolved), 7 reactions (Lund + 2).
        @test length(m.components) == 9
        @test [c.name for c in particles(m)] == ["HET", "PHO", "POM", "PAT"]
        @test [c.name for c in liquids(m)] == ["O2", "IC", "NH4", "HPO4", "DOM"]
        @test [r.name for r in m.reactions][6:7] == ["Inactivation", "Bacterivory"]

        # thesis_model globals/units differ from modelLund (τ = 1e-3, ρ ≈ 1).
        @test m.osmosis_rate ≈ 1.0e-3
        @test m.cohesion_submodel.zeta_0 ≈ 1.0
        @test particles(m)[1].density ≈ 1.2
        @test particles(m)[4].transport_rate ≈ 1.0e-6   # PAT
        @test particles(m)[4].attachment_sand ≈ 2.0

        # Inactivation: first-order PAT die-off, inert in the flowing phase.
        inact = m.reactions[6]
        @test inact.order == Dict("PAT" => 1.0)
        @test inact.stoichiometric_coefficients == Dict("PAT" => -1.0)
        @test inact.efficiency_flowing == 0.0
        @test inact.nominal_rate ≈ 0.4

        # Bacterivory: first-order PAT, HET-Monod (kPred = 0.002), returns HET;
        # the only flowing-active reaction, scaled by water_factor.
        bac = m.reactions[7]
        @test bac.order == Dict("PAT" => 1.0)
        @test bac.half_saturation_constants == Dict("HET" => 0.002)   # Monod on HET, not DOM
        @test bac.stoichiometric_coefficients == Dict("PAT" => -1.0, "HET" => 1.0)
        @test bac.nominal_rate ≈ 20.0
        @test bac.efficiency_flowing ≈ 1.0e-3                          # default water_factor

        # All Lund reactions are inert in the flowing phase (r_water = 0).
        @test all(m.reactions[i].efficiency_flowing == 0.0 for i in 1:5)

        # Constructor knobs propagate.
        m2 = modelPathogen(water_factor=5e-3, sand_pathogen=0.1, dark_respiration=0.05)
        @test m2.reactions[7].efficiency_flowing ≈ 5e-3
        @test particles(m2)[4].sand_attachment_factor ≈ 0.1
        @test m2.reactions[2].minimum_light_factor ≈ 0.05             # phototroph dark resp.
        # Default PAT does not preferentially attach to sand.
        @test particles(m)[4].sand_attachment_factor == 0.0

        # The generic reaction kernel accepts the model (one order term each).
        @test MPCSSF._reaction_kernel(m) isa NamedTuple
    end

    @testset "pathogen reaction kinetics" begin
        # Directly evaluate the two new reactions through the generic kernel, at
        # 20 °C where μ reduces to the nominal rates, to confirm the encoding
        # reproduces run_pathogen's Inactivation (μ·PAT) and Bacterivory
        # (μ·PAT·HET/(HET+kPred)) terms.
        m = modelPathogen()
        comps = vcat(particles(m), liquids(m))
        kernel = MPCSSF._reaction_kernel(m)
        mu = compute_reaction_rates(m, 20)                 # μ₂₀ = nominal at 20 °C
        @test mu[6] ≈ 0.4 && mu[7] ≈ 20.0

        HET, PAT = 1.0, 2.0
        X = [HET 0.0 0.5 PAT]                              # HET PHO POM PAT
        S = [0.1 0.1 0.1 0.1 0.1]                          # O2 IC NH4 HPO4 DOM
        local_ = MPCSSF._local(X, S, kernel.num_idx, kernel.den_idx)
        light = ones(1, length(mu))
        rx = MPCSSF._evaluate_reactions(local_, kernel, [1.0], mu, light)

        kPred = 0.002
        @test rx[1, 6] ≈ 0.4 * PAT                                    # inactivation
        @test rx[1, 7] ≈ 20.0 * PAT * HET / (HET + kPred) rtol=1e-9   # bacterivory
    end

    @testset "pathogen simulate (smoke)" begin
        # The 4-particle / 7-reaction pathogen model runs end-to-end through the
        # generalized solver (efficiency_flowing + sand_attachment_factor paths).
        f = addgridpoints(SandFilter(), 20)
        m = modelPathogen()
        inflow = Float64[1e-3, 1e-3, 0.0, 1e-5, 1e-2, 1e-2, 1e-5, 0.0, 1e-4] ./ 10
        r = simulate(State(f, m); inflow_concentrations=inflow,
                     simulation_time=1e-6, time_step=1e-8, n_frames=5, quiet=true)
        N = length(f.grid.centers)
        @test r.flag == "OK"
        @test size(r.frames[:concentration_biofilm]) == (N, 5, 13)    # 2kP+kL
        @test size(r.frames[:concentration_flowing]) == (N, 5, 9)     # kP+kL
        @test reaction_names(r)[6:7] == ["Inactivation", "Bacterivory"]
        @test size(reaction_rates(r), 3) == 7
    end

    @testset "light-factor floor (dark respiration)" begin
        # The light factor is a FLOOR: max(fdark, I_eff·e^{1-I_eff}), the port of
        # the authoritative slow-sand-filtration form. This is the reconciliation
        # away from the old additive max(0, fdark + ...) that MPC-SSF/older Julia
        # used. The cross-implementation float floor in a full sim masks this
        # difference, so validate the formula directly.
        le = [0.0, 0.05, 0.2, 0.8]            # effective light across depth
        fdark = 0.1
        lf = MPCSSF._light_factor_floor(le, [fdark])
        @test vec(lf) ≈ max.(fdark, le)                  # floor form
        @test vec(lf) == [0.1, 0.1, 0.2, 0.8]            # both branches exercised
        @test vec(lf) != fdark .+ le                     # NOT the additive form
        # fdark = 0 reduces to the bare light term (and matches the old additive).
        @test vec(MPCSSF._light_factor_floor(le, [0.0])) ≈ le
        # one column per light-dependent reaction
        lf2 = MPCSSF._light_factor_floor(le, [0.1, 0.3])
        @test size(lf2) == (4, 2)
        @test lf2[:, 1] ≈ max.(0.1, le) && lf2[:, 2] ≈ max.(0.3, le)
        # reproduces the full authoritative chain: I_eff·e^{1-I_eff} then floor
        ieff = 0.8 .* exp.(-[0.0, 0.5, 2.0]) ./ 1.08
        le3 = ieff .* exp.(1 .- ieff)
        @test vec(MPCSSF._light_factor_floor(le3, [0.1])) ≈ max.(0.1, le3)
    end

    @testset "adaptive CFL time-stepping" begin
        f = addgridpoints(SandFilter(), 20)
        m = modelLund()
        inflow = Float64[1e-2, 1e-3, 0.0, 1e-4, 1e-2, 5e-3, 4e-3, 1e-4, 1.0]

        # Adaptive run completes; dt ramps up so many steps fit a long horizon
        # that a fixed 1e-10 step never could in the same wall time.
        r = simulate(State(f, m); inflow_concentrations=inflow, simulation_time=1e-5,
                     time_step=:adaptive, adaptive_initial_dt=1e-10,
                     adaptive_max_dt=1e-6, cfl_factor=0.99, n_frames=5, quiet=true)
        @test r.flag == "OK"
        @test r.simulation_data[:time_step] == "adaptive"

        # Capping the adaptive step at AdaptiveMaxDt reproduces the fixed step
        # exactly (dt is clamped to the cap every iteration).
        rF = simulate(State(f, m); inflow_concentrations=inflow, simulation_time=3e-10,
                      time_step=1e-10, n_frames=4, quiet=true)
        rA = simulate(State(f, m); inflow_concentrations=inflow, simulation_time=3e-10,
                      time_step=:adaptive, adaptive_initial_dt=1e-10,
                      adaptive_max_dt=1e-10, n_frames=4, quiet=true)
        @test rF.flag == rA.flag
        @test times(rF) == times(rA)
        @test concentration(rF, "HET", :flowing) == concentration(rA, "HET", :flowing)
        @test concentration(rF, "DOM", :flowing) == concentration(rA, "DOM", :flowing)
    end

    @testset "run chaining (resume + concatenate)" begin
        f = addgridpoints(SandFilter(), 20)
        m = simpleModel()
        inflow = Float64[1e-2, 0.0, 1.0]
        T = 5e-4; dt = 1e-5   # T is an exact multiple of dt so the grids align

        rfull = simulate(State(f, m); inflow_concentrations=inflow,
                         simulation_time=2T, time_step=dt, n_frames=5, quiet=true)
        r1 = simulate(State(f, m); inflow_concentrations=inflow,
                      simulation_time=T, time_step=dt, n_frames=5, quiet=true)
        s = final_state(r1)
        @test s.time ≈ r1.time_final
        @test s isa State
        r2 = simulate(s; inflow_concentrations=inflow, simulation_time=T,
                      time_step=dt, n_frames=5, quiet=true)

        # resuming reproduces the continuous run exactly (T is on the dt grid)
        @test concentration(rfull, "Microorganism", :flowing)[:, end] ==
              concentration(r2, "Microorganism", :flowing)[:, end]
        @test concentration(rfull, "Nutrient", :flowing)[:, end] ==
              concentration(r2, "Nutrient", :flowing)[:, end]

        # concatenate joins the two runs, dropping r2's duplicated first frame
        rcat = r1 + r2
        @test rcat.time_start == r1.time_start
        @test rcat.time_final == r2.time_final
        @test rcat.flag == "OK"
        @test length(times(rcat)) == length(times(r1)) + length(times(r2)) - 1
        @test issorted(times(rcat)) && allunique(times(rcat))
        @test size(concentration(rcat, "Microorganism", :flowing), 2) == length(times(rcat))

        # incompatible joins error
        @test_throws Exception concatenate(r1, r1)   # r1.time_start ≠ r1.time_final
        # UNINITIATED is the identity element
        blank = Results(f, m)
        @test (blank + r1).flag == r1.flag
        @test (r1 + blank).flag == r1.flag
    end

    # Golden-master parity vs MATLAB. Runs against the committed reference
    # (test/golden/reference/, produced by export_reference.m) by default; set
    # MPCSSF_GOLDEN_REF to compare against a freshly exported reference instead.
    # See test/golden/README.md.
    @testset "golden-master vs MATLAB (fixed step)" begin
        include(joinpath(@__DIR__, "golden", "compare.jl"))
        refdir = get(ENV, "MPCSSF_GOLDEN_REF", GOLDEN_DEFAULT_REF)
        r = golden_compare(refdir; verbose=true)
        @test r.flags_agree
        @test r.match
    end

    # Adaptive-CFL golden-master: modelLund with TimeStep=:adaptive vs the MATLAB
    # reference (test/golden/reference_adaptive/, from export_adaptive_reference.m).
    @testset "golden-master vs MATLAB (adaptive)" begin
        include(joinpath(@__DIR__, "golden", "compare_adaptive.jl"))
        refdir = get(ENV, "MPCSSF_GOLDEN_ADAPTIVE_REF", GOLDEN_ADAPTIVE_REF)
        r = golden_compare_adaptive(refdir; verbose=true)
        @test r.flags_agree
        @test r.match
        # the adaptive stepping is deterministic here: identical step count
        @test r.jl_steps == r.ref_steps
    end

    # Pathogen golden-master: modelPathogen vs the AUTHORITATIVE slow-sand-
    # filtration @SDfilter/run_pathogen.m + thesis_model.mat (reference exported
    # by export_pathogen_reference.m into test/golden/reference_pathogen/). A
    # seeded mature biofilm makes the biofilm-phase reactions (death, hydrolysis,
    # pathogen inactivation, bacterivory) fire; the run also exercises
    # sand_pathogen (differential PAT→sand attachment) and the flowing-phase
    # water_factor. See test/golden/README.md.
    @testset "golden-master vs MATLAB (pathogen)" begin
        include(joinpath(@__DIR__, "golden", "compare_pathogen.jl"))
        refdir = get(ENV, "MPCSSF_GOLDEN_PATHOGEN_REF", PGOLD_DEFAULT_REF)
        r = golden_compare_pathogen(refdir; verbose=true)
        @test r.flags_agree
        @test r.match
    end

    # Light-active pathogen golden: constant light + dark_respiration>0 with a
    # phototroph-heavy seed, vs run_pathogen.m. End-to-end check that light-active
    # runs reproduce the authoritative code (the dark-respiration light floor is
    # unit-tested separately — the cross-implementation float floor here, ~3e-7,
    # is reaction-driven and masks the light-form difference). atol relaxed to
    # 1e-6 accordingly. See test/golden/README.md.
    @testset "golden-master vs MATLAB (pathogen, light+dark)" begin
        refdir = get(ENV, "MPCSSF_GOLDEN_PATHOGEN_LIGHT_REF", PGOLD_LIGHT_REF)
        r = golden_compare_pathogen(refdir; run_kwargs=PGOLD_LIGHT_RUN,
                                    atol=1e-6, verbose=true)
        @test r.flags_agree
        @test r.match
    end

    # None of the golden suites above can see a change in the cohesion potential:
    # two start from phi_b = 0, where every form vanishes, and the pathogen pair
    # run 1e-5 days. These check the wiring the goldens cannot.
    include(joinpath(@__DIR__, "cohesion_discrimination.jl"))
end
