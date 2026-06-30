# Port of src/@SandFilter/{SandFilter,addGridPoints,computePorosity,
#                            getCahnHilliardMatrices}.m
#
# Filter geometry, 1D staggered grid, porosity profile, light attenuation, and
# Cahn-Hilliard matrix assembly. MATLAB stored the grid as a struct
# (Boundaries + Centers) with Dependent properties GridSize/GridZero and the
# light-attenuation eta vectors; here those are a `Grid` struct plus functions.

"Staggered 1D grid: cell boundaries (faces) and centers (port of GridPoints)."
struct Grid
    boundaries::Vector{Float64}
    centers::Vector{Float64}
end

"Build a `Grid` from face positions; centers are face midpoints (set.GridPoints)."
function Grid(boundaries::AbstractVector{<:Real})
    b = collect(float.(boundaries))
    centers = 0.5 .* (b[1:end-1] .+ b[2:end])
    return Grid(b, centers)
end

"Default diel light forcing (port of the anonymous fn in SandFilter.m)."
default_light_irradiation(t) = 0.8 * max(sin(2π * (t - 13 / 48)) + 31 / 50, 0) / (1 + 31 / 50)

"""
    SandFilter(; height=1.0, depth=1.0, sand_porosity=0.4, sand_roughness=5e-3,
                 inflow_velocity=0.3*24, temperature=15+273,
                 light_irradiation=default_light_irradiation,
                 light_attenuation_water=0.32, light_attenuation_sand=1500.0,
                 grid=nothing)

Filter geometry and discretization. Defaults mirror `SandFilter.m`. `depth` is
the packed-sand depth (z ≥ 0) and `height` the supernatant water height
(z ≤ 0); the sand surface is at z = 0.
"""
Base.@kwdef mutable struct SandFilter
    height::Float64 = 1.0
    depth::Float64 = 1.0
    sand_porosity::Float64 = 0.4
    sand_roughness::Float64 = 5e-3
    inflow_velocity::Float64 = 0.3 * 24
    temperature::Float64 = 15 + 273
    light_irradiation::Function = default_light_irradiation
    light_attenuation_water::Float64 = 0.32
    light_attenuation_sand::Float64 = 1500.0
    grid::Union{Grid,Nothing} = nothing
end

# --- grid accessors ----------------------------------------------------------

# mean(diff(boundaries)) telescopes to (last - first)/(n - 1) for any grid.
"Mean cell size dz (port of get.GridSize)."
function gridsize(f::SandFilter)
    b = f.grid.boundaries
    return (b[end] - b[1]) / (length(b) - 1)
end

"Index of the cell straddling the sand surface z = 0 (port of get.GridZero)."
function gridzero(f::SandFilter)
    dz = gridsize(f)
    i = findfirst(z -> abs(z) < dz / 2, f.grid.centers)
    i === nothing && error("no cell center within dz/2 of z = 0; check the grid")
    return i
end

"""
    addgridpoints(f::SandFilter, n::Integer) -> SandFilter

Port of `addGridPoints.m`: build the staggered grid with `n` interior cells.
`dz = depth/(n + 1/2)`; faces run from −height up to depth + dz. If a face lands
on z ≈ 0 the cell count is bumped by one (so no computational boundary sits on
the sand surface). Sets `f.grid` and returns `f`.
"""
function addgridpoints(f::SandFilter, n::Integer)
    dz = f.depth / (n + 1 / 2)
    decreasing = collect(f.depth:-dz:-f.height)      # depth → ≥ −height
    lower = vcat(f.depth + dz, decreasing)
    if all(abs.(lower) .> 1e-16)
        f.grid = Grid(reverse(lower))                # ascending faces
        return f
    else                                             # a face ≈ z = 0; retry
        return addgridpoints(f, n + 1)
    end
end

# --- porosity ----------------------------------------------------------------

"""
    computeporosity(f::SandFilter, z) -> porosity

Port of `computePorosity.m`. Linear ramp from water (porosity → 1 for z < 0)
to packed sand (porosity → `sand_porosity` for z ≥ 0), clamped to
`[sand_porosity, 1]`. `z` may be a scalar or a vector.
"""
function computeporosity(f::SandFilter, z)
    ε0 = f.sand_porosity
    δ = f.sand_roughness
    return clamp.((ε0 - 1) / δ .* z .+ ε0, ε0, 1.0)
end

# --- light attenuation (port of the Dependent eta properties) ---------------

"Optical depth η in the supernatant water (port of LightAttenuationEtaWater)."
function light_attenuation_eta_water(f::SandFilter)
    centers = f.grid.centers
    return f.light_attenuation_water .* (centers .+ f.height)
end

"Optical depth η through the sand roughness layer and packed bed (port of LightAttenuationEtaSand)."
function light_attenuation_eta_sand(f::SandFilter)
    centers = f.grid.centers
    δ = f.sand_roughness
    ε0 = f.sand_porosity
    eta = zeros(length(centers))
    rough = (centers .>= -δ) .& (centers .< 0)       # roughness layer [−δ, 0)
    packed = centers .>= 0                            # packed bed z ≥ 0
    cr = centers[rough]
    eta[rough] = (1 - ε0) .* (cr .+ δ) .* (1 .+ (cr .- δ) ./ (2δ))
    eta[packed] = (1 - ε0) .* (δ / 2 .+ centers[packed])
    return f.light_attenuation_sand .* eta
end

# --- Cahn-Hilliard matrix assembly ------------------------------------------

"""
    Triplets(I, J, V, m, n)

COO sparse-matrix triplets. `SparseArrays.sparse(t)` materializes the `m × n`
matrix (duplicate (I,J) entries are summed, matching MATLAB `sparse`).
"""
struct Triplets
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{Float64}
    m::Int
    n::Int
end
SparseArrays.sparse(t::Triplets) = sparse(t.I, t.J, t.V, t.m, t.n)

"""
    get_cahn_hilliard_matrices(f::SandFilter, model, upwinded::Bool=false)
        -> (convection, diffusion, diffusion_mobility)

Port of `getCahnHilliardMatrices.m`. Returns three [`Triplets`](@ref) (each
`2n0 × 2n0`, where `n0 = gridzero(f)`) for the convection, diffusion, and
mobility-weighted diffusion operators of the Cahn-Hilliard system ("SOLVER A").
`model` must expose `model.cohesion_submodel.kappa`.
"""
function get_cahn_hilliard_matrices(f::SandFilter, model, upwinded::Bool=false)
    centers = f.grid.centers
    boundaries = f.grid.boundaries
    porosity_centers = computeporosity(f, centers)
    porosity_boundaries = computeporosity(f, boundaries)

    dz = gridsize(f)
    n0 = gridzero(f)
    N = 2n0

    # 2×2 local stencil flattened column-major: rows [1,2,1,2], cols [1,1,2,2].
    lr = (1, 2, 1, 2)
    lc = (1, 1, 2, 2)
    conv_stencil = upwinded ? (-1.0, 1.0, 0.0, 0.0) : (-0.5, 0.5, -0.5, 0.5)
    diff_stencil = (1.0, -1.0, -1.0, 1.0)

    q = f.inflow_velocity
    κ = model.cohesion_submodel.kappa

    # CONVECTION: block rows/cols 1..n0, n0-1 cell pairs.
    cI = Int[]; cJ = Int[]; cV = Float64[]
    for k in 1:(n0 - 1)               # offset = k - 1
        for a in 1:4
            r = lr[a] + (k - 1)
            push!(cI, r)
            push!(cJ, lc[a] + (k - 1))
            push!(cV, q * conv_stencil[a] / porosity_centers[r] / dz)
        end
    end
    convection = Triplets(cI, cJ, cV, N, N)

    # DIFFUSION: clamped indices, n0+1 columns.
    dI = Int[]; dJ = Int[]; dV = Float64[]
    for k in 1:(n0 + 1)               # offset value = k - 1
        for a in 1:4
            r = clamp(lr[a] + (k - 1) - 1, 1, n0)
            c = clamp(lc[a] + (k - 1) - 1, 1, n0)
            push!(dI, r)
            push!(dJ, c)
            push!(dV, κ * diff_stencil[a] / dz^2)
        end
    end
    diffusion = Triplets(dI, dJ, dV, N, N)

    # MOBILITY-WEIGHTED DIFFUSION: columns shifted into the second block (+n0);
    # weight λ = porosity at interior boundaries 2..n0.
    mI = Int[]; mJ = Int[]; mV = Float64[]
    for k in 1:(n0 - 1)               # offset = k - 1
        λ = porosity_boundaries[k + 1]
        for a in 1:4
            r = lr[a] + (k - 1)
            push!(mI, r)
            push!(mJ, lc[a] + (k - 1) + n0)
            push!(mV, -(diff_stencil[a] * λ) / porosity_centers[r] / dz^2)
        end
    end
    diffusion_mobility = Triplets(mI, mJ, mV, N, N)

    return (convection=convection, diffusion=diffusion,
            diffusion_mobility=diffusion_mobility)
end
