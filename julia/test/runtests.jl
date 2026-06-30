using MPCSSF
using Test

# Phase 1 verification scaffold. As each unit is ported, replace the
# `@test_broken`/`error` placeholders with real assertions, ideally golden-master
# comparisons against MATLAB reference outputs exported for SimpleModel and
# LundMPCModel (see ../.claude/julia-port-outline.md §5).

@testset "MPCSSF" begin
    @testset "construction" begin
        # Smoke tests: structs build with sensible defaults.
        p = Particle(name="HET", density=1.0)
        @test p.name == "HET"
        @test p.density == 1.0

        l = Liquid(name="DOM", density=1000.0)
        @test l isa Component

        f = SandFilter()
        @test f.sand_porosity == 0.4
        @test f.light_irradiation(0.0) isa Real

        m = Model(components=[p, l])
        @test length(MPCSSF.particles(m)) == 1
        @test length(MPCSSF.liquids(m)) == 1
    end

    @testset "ecological layer (to port)" begin
        rx = Reaction(nominal_rate=1.0)
        @test_throws Exception MPCSSF.compute_rate(rx, 293)
    end

    @testset "solver (to port)" begin
        # Will become a real short-run golden-master test.
        f = SandFilter()
        m = Model()
        @test_throws Exception State(f, m)
    end
end
