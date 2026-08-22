# Port of src/@State/simulate.m  (the core fixed-step solver)
#
# Ported block-by-block, preserving the MATLAB section structure so it can be
# diffed against the source. Phase 1 implements the FIXED-step path only (parity
# with current src/). Adaptive (CFL) time-stepping and the pathogen variant are
# later phases.
#
# Sections, mirroring simulate.m:
#   I/II.  parameters (filter geometry, model params, reaction kernel)
#   III.   initial conditions (volume fractions, initial velocities)
#   IV/V.  pre-allocation + output scaffolding
#   VI.    time loop:
#            - global/local concentrations & volume fractions
#            - clogging guard
#            - reaction terms (light, attachment/detachment, transfer, eco)
#            - SOLVER A: Cahn-Hilliard implicit solve -> biofilm velocity
#            - SOLVER B: finite-volume update (upwind convection + dispersion)
#            - negativity guards, frame capture

const _REALMIN = floatmin(Float64)   # MATLAB realmin

# --- reaction kernel (port of the listK / listOrder setup) ------------------
# Precomputes, per reaction: the columns of the local-concentration matrix that
# carry a half-saturation (Monod) term and their K values, plus the column of
# the order-1 species. The local matrix is [X (particles) | S (liquids) |
# quotient columns]; this requires components ordered particles-then-liquids
# (true for the Lund/Rosenqvist presets), so we build over that order.
function _reaction_kernel(model::Model)
    comps = vcat(particles(model), liquids(model))   # particles-then-liquids
    rxs = model.reactions
    nRx = length(rxs)
    hsK = lookup_half_saturation_constants(rxs, comps)   # nComp × nRx (NaN absent)
    ord = lookup_order(rxs, comps)                       # nComp × nRx
    quot = lookup_quotients(rxs, comps)                  # K nQuot×nRx, num/den idx
    nComp = length(comps)
    listK = vcat(hsK, quot.K)                            # (nComp+nQuot) × nRx

    monod_cols = [findall(!isnan, @view listK[:, i]) for i in 1:nRx]
    monod_K = [listK[monod_cols[i], i] for i in 1:nRx]

    order_col = Vector{Int}(undef, nRx)
    for i in 1:nRx
        r = findfirst(!=(0), @view ord[:, i])
        r === nothing &&
            error("reaction $i has no kinetic-order term; simulate requires exactly one per reaction")
        order_col[i] = r
    end
    return (; monod_cols, monod_K, order_col,
            num_idx=quot.num_idx, den_idx=quot.den_idx)
end

# Build a region's local-concentration matrix [X | S | quotients].
function _local(X, S, num_idx, den_idx)
    XS = hcat(X, S)
    isempty(num_idx) && return XS
    Q = XS[:, num_idx] ./ (XS[:, den_idx] .+ _REALMIN)
    return hcat(XS, Q)
end

# Dark-respiration light floor. Port of the authoritative slow-sand-filtration
# form `I = max(fdark, I_eff·e^{1-I_eff})` (@SDfilter/run_biofilm.m,
# run_pathogen.m): the light factor is at least `min_light` even in darkness — a
# FLOOR, not an additive offset. The per-reaction `minimum_light_factor` stands
# in for the global `dark_respiration` (only light-dependent reactions carry it).
#
# `light_effective` is the per-cell I_eff·e^{1-I_eff} (length nCells); `min_light`
# is the floor per light-dependent reaction (length nLightDep). Returns an
# nCells × nLightDep matrix. NOTE: MPC-SSF's simulate.m and older SSF.jl used
# the *additive* form `max(0, min_light + light_effective)`, which double-counts
# the baseline at high light; this floor form supersedes it (see test/golden/README.md).
_light_factor_floor(light_effective::AbstractVector, min_light::AbstractVector) =
    max.(reshape(min_light, 1, :), light_effective)

# Adaptive-CFL source bounds (Diehl2025 Theorem 1 / Section 3.7). These are the
# per-step reaction contributions to the time-step bound. They are typed
# top-level functions rather than closures defined inside the time loop: as
# closures capturing the loop's variables they were boxed/type-unstable and
# became the single largest per-step hotspot. The numerics are unchanged
# (identical loop order and arithmetic), so golden masters stay bit-identical.

# Liquid-consumption bound for regions Le/Lf: max over cells n and liquids l of
# |Σ_{j∈{growth}} σ_L[l,j]·μ_j·X[n,j] / (S[n,l] + K[l,j])|. Only the two growth
# reactions (particle columns 1, 2) enter the bound.
function _cfl_liquid_bound(X::AbstractMatrix, S::AbstractMatrix,
                           sigmaL::AbstractMatrix, muRates::AbstractVector,
                           K_CFL::AbstractMatrix, kL::Integer)
    m = 0.0
    @inbounds for n in axes(X, 1), l in 1:kL
        acc = 0.0
        for j in 1:2
            acc += sigmaL[l, j] * muRates[j] * X[n, j] / (S[n, l] + K_CFL[l, j])
        end
        m = max(m, abs(acc))
    end
    return m
end

# Hydrolysis quotient bound for the particle regions: max over cells of
# X_POM / (X_POM + K_hyd·X_HET) (columns 3 and 1); 0/0 ⇒ 0 (matches NaN→0).
function _cfl_hydrolysis_bound(X::AbstractMatrix, K_Hyd::Real)
    m = 0.0
    @inbounds for n in axes(X, 1)
        den = X[n, 3] + K_Hyd * X[n, 1]
        xi = den == 0 ? 0.0 : X[n, 3] / den
        m = max(m, xi)
    end
    return m
end

# In-place port of the evaluateReactions subfunction. Writes the nCells×nRx
# reaction-rate matrix into `rx` (preallocated in the time loop). The scalar loop
# computes the Monod minimum without per-reaction temporaries; the arithmetic and
# reduction order are identical to the broadcast form, so results are bit-for-bit
# unchanged (monod terms are ≤ 1, so seeding the min with 1.0 is a no-op).
function _evaluate_reactions!(rx, local_, kernel, phi, muRates, lightFactor)
    nCells, nRx = size(rx)
    @inbounds for i in 1:nRx
        cols = kernel.monod_cols[i]
        Ki = kernel.monod_K[i]
        oc = kernel.order_col[i]
        mu_i = muRates[i]
        for n in 1:nCells
            monod = 1.0
            for kk in eachindex(cols)
                x = local_[n, cols[kk]]
                monod = min(monod, (x + _REALMIN) / (Ki[kk] + x + _REALMIN))
            end
            rx[n, i] = phi[n] * mu_i * lightFactor[n, i] * monod * local_[n, oc]
        end
    end
    return rx
end

# Allocating wrapper (external callers / tests). `phi` is per-cell; returns nCells×nRx.
_evaluate_reactions(local_, kernel, phi, muRates, lightFactor) =
    _evaluate_reactions!(zeros(size(local_, 1), length(muRates)),
                         local_, kernel, phi, muRates, lightFactor)

_meanval(x) = sum(x) / length(x)

"""
    simulate(state::State; inflow_concentrations=nothing, simulation_time=1.0,
             time_step=1e-5, n_frames=200, clogging_fraction=0.99,
             upwinded=false, quiet=false,
             cfl_factor=0.99, adaptive_velocity_factor=0.0,
             adaptive_time_tolerance=1e-3, adaptive_initial_dt=5e-7,
             adaptive_max_dt=5e-7) -> Results

Run the slow-sand-filtration simulation from `state`. `time_step` is either a
number (fixed step) or `:adaptive`, which recomputes the step each iteration
from a per-step CFL bound (ported from the MATLAB `simulate` adaptive path,
itself from slow-sand-filtration `@SDfilter/run_biofilm.m`). The adaptive path
grows `dt` by at most `(1 + adaptive_time_tolerance)` per step toward
`cfl_factor/max(w)`, capped at `adaptive_max_dt`, starting from
`adaptive_initial_dt`; it assumes the Lund/Rosenqvist model structure (see
[`modelLund`](@ref)). `inflow_concentrations` may be `nothing` (zero inflow), a
vector of length `kP+kL`, or a function `t -> vector`. Returns a
[`Results`](@ref) whose `flag` is `"OK"`, `"CLOGGED"`, `"BIOFILM"`, or
`"FLOWING"`.
"""
function simulate(state::State;
                  inflow_concentrations=nothing,
                  simulation_time::Real=1.0,
                  time_step=1e-5,
                  n_frames::Integer=200,
                  clogging_fraction::Real=0.99,
                  upwinded::Bool=false,
                  quiet::Bool=false,
                  cfl_factor::Real=0.99,
                  adaptive_velocity_factor::Real=0.0,
                  adaptive_time_tolerance::Real=1e-3,
                  adaptive_initial_dt::Real=5e-7,
                  adaptive_max_dt::Real=5e-7,
                  implicit_osmosis::Bool=false)
    f = state.filter
    model = state.model
    temperature = f.temperature

    # ---- I. filter geometry ------------------------------------------------
    centers = f.grid.centers
    boundaries = f.grid.boundaries
    N = length(centers)
    dz = gridsize(f)
    n0 = gridzero(f)
    porosity_centers = computeporosity(f, centers)         # length N
    porosity_boundaries = computeporosity(f, boundaries)   # length N+1

    mats = get_cahn_hilliard_matrices(f, model, upwinded)
    Smat = sparse(mats.convection)
    DDmat = sparse(mats.diffusion)
    CH0 = spdiagm(0 => ones(2n0)) - DDmat
    Dmob = mats.diffusion_mobility                         # triplets, scaled per step

    eta_water = light_attenuation_eta_water(f)
    eta_sand = light_attenuation_eta_sand(f)
    volumetric_flow = f.inflow_velocity
    volume_avg_velocity = volumetric_flow ./ porosity_boundaries   # length N+1

    # ---- II. model parameters ---------------------------------------------
    P = particles(model)
    L = liquids(model)
    kP = length(P)
    kL = length(L)
    densityP = _meanval([p.density for p in P])            # assumption 1
    densityL = _meanval([l.density for l in L])
    attachment_rates = reshape([p.attachment_sand for p in P], 1, kP)
    transport_particle_rates = reshape([p.transport_rate for p in P], 1, kP)
    transport_liquid_rates = reshape([l.transport_rate for l in L], 1, kL)
    alpha = reshape(vcat([p.dispersivity for p in P], [l.dispersivity for l in L]), 1, kP + kL)
    beta = model.biofilm_porosity
    tau = model.osmosis_rate
    sigmaP = stoichiometric_matrix_particles(model)        # kP × nRx
    sigmaL = stoichiometric_matrix_liquids(model)          # kL × nRx
    dpsi_fun = model.cohesion_submodel.potential_gradient
    zeta_0 = model.cohesion_submodel.zeta_0
    mobility = model.cohesion_submodel.mobility

    muRates = compute_reaction_rates(model, temperature)   # length nRx
    nRx = length(muRates)
    # Per-reaction phase efficiencies (default 1.0). `efficiency_flowing` is the
    # port of the MATLAB `r_water` scaling: it multiplies each reaction's rate in
    # the flowing phase (the pathogen model zeroes all but bacterivory and scales
    # that by `water_factor`). `efficiency_biofilm` scales the biofilm/enclosed
    # phases. Both default to 1.0, leaving non-pathogen models unchanged.
    eff_biofilm = reshape(Float64[r.efficiency_biofilm for r in model.reactions], 1, nRx)
    eff_flowing = reshape(Float64[r.efficiency_flowing for r in model.reactions], 1, nRx)
    # Per-particle sand-attachment scaling (default 1.0); see `Particle`.
    sand_factors = reshape(Float64[p.sand_attachment_factor for p in P], 1, kP)
    light_optimal = isempty(model.reactions) ? 1.0 :
                    maximum(r.optimal_light_factor for r in model.reactions)
    attenuation_particles = reshape([p.attenuation for p in P], 1, kP)
    light_dep = [r.is_light_dependent for r in model.reactions]
    min_light = [r.minimum_light_factor for r in model.reactions][light_dep]
    light_inh = [r.light_inhibition for r in model.reactions]
    inh_dep = light_inh .> 0
    comp_dep = [r.is_light_complement for r in model.reactions]
    # Inhibited-only models have no optimal_light_factor; fall back to 1.0 so the
    # intensity normalization below stays finite (value unused otherwise).
    light_optimal = light_optimal > 0 ? light_optimal : 1.0
    kernel = _reaction_kernel(model)

    # ---- III. initial conditions ------------------------------------------
    t0 = state.time
    globalBiofilm = global_concentration_biofilm(state)    # N × (2kP+kL)
    globalFlowing = global_concentration_flowing(state)    # N × (kP+kL)
    phiW = copy(state.enclosed_water_volume)               # N
    velBiofilm = copy(state.velocity.biofilm)              # N-1

    phiBiofilm = phiW .+ vec(sum(@view(globalBiofilm[:, 1:kP]), dims=2)) ./ densityP .+
                 vec(sum(@view(globalBiofilm[:, kP+1:2kP]), dims=2)) ./ densityP .+
                 vec(sum(@view(globalBiofilm[:, 2kP+1:2kP+kL]), dims=2)) ./ densityL
    phiBdy = 0.5 .* (phiBiofilm[2:end] .+ phiBiofilm[1:end-1])
    velFlowing = vcat(volume_avg_velocity[1],
                      (volume_avg_velocity[2:end-1] .- velBiofilm .* phiBdy) ./ (1 .- phiBdy),
                      volume_avg_velocity[end])            # length N+1

    # ---- IV/V. pre-allocation + output ------------------------------------
    numFrames = n_frames
    timeFrames = zeros(numFrames)
    concFramesBiofilm = zeros(N, numFrames, 2kP + kL)
    concFramesFlowing = zeros(N, numFrames, kP + kL)
    concFramesWater = zeros(N, numFrames)
    velFramesBiofilm = zeros(N - 1, numFrames)
    velFramesFlowing = zeros(N - 1, numFrames)
    concFramesReaction = zeros(N, numFrames, nRx)   # biofilm eco reaction rates

    timeSnap = collect(range(t0, t0 + simulation_time; length=numFrames))
    counter = 1
    timeFrames[counter] = t0
    concFramesBiofilm[:, counter, :] = globalBiofilm
    concFramesFlowing[:, counter, :] = globalFlowing
    concFramesWater[:, counter] = densityL .* phiW
    velFramesBiofilm[:, counter] = velBiofilm
    velFramesFlowing[:, counter] = velFlowing[2:end-1]
    counter += 1

    results = Results(f, model)
    results.time_start = t0
    results.flag = "OK"

    # inflow handling
    inflow_fn = if inflow_concentrations === nothing
        _t -> zeros(1, kP + kL)
    elseif inflow_concentrations isa Function
        _t -> reshape(collect(float.(inflow_concentrations(_t))), 1, kP + kL)
    else
        row = reshape(collect(float.(inflow_concentrations)), 1, kP + kL)
        _t -> row
    end

    # ---- VI. time integration ---------------------------------------------
    # Fixed step (numeric time_step) or per-step CFL bound (time_step=:adaptive).
    is_adaptive = time_step === :adaptive
    dt = is_adaptive ? float(adaptive_initial_dt) : float(time_step)
    if is_adaptive
        # CFL constants (Lund-structured model assumption; see docstring).
        alphaP = alpha[1]
        alphaL = alpha[kP+1]
        hs = half_saturation_constants(model)              # nComp × nRx, NaN absent
        K_HetGrowth = [isnan(x) ? Inf : x for x in hs[kP+1:kP+kL, 1]]
        K_PhoGrowth = [isnan(x) ? Inf : x for x in hs[kP+1:kP+kL, 2]]
        K_CFL = hcat(K_HetGrowth, K_PhoGrowth)             # kL × 2
        qvals = [x for x in quotients(model).K if !isnan(x)]   # hydrolysis POM/HET
        K_Hyd = isempty(qvals) ? Inf : maximum(qvals)
    end
    results.simulation_data[:time_step] = is_adaptive ? "adaptive" : dt
    t = t0
    step_times = Float64[t0]        # t after each step; diff gives the dt sequence
    # reused per-step reaction buffers (overwritten each step by _evaluate_reactions!)
    ecoBiofilmBuf  = zeros(N, nRx)
    ecoEnclosedBuf = zeros(N, nRx)
    ecoFlowingBuf  = zeros(N, nRx)
    # light factor: non-light-dependent columns are 1 for all time; only the
    # light-dependent columns are rewritten each step, so allocate/fill once.
    lightFactor = ones(N, nRx)
    while t < t0 + simulation_time
        globalConcInflow = inflow_fn(t)

        # global concentrations by region
        globalMatrix    = globalBiofilm[:, 1:kP]
        globalEnclosedP = globalBiofilm[:, kP+1:2kP]
        globalEnclosedL = globalBiofilm[:, 2kP+1:2kP+kL]
        globalFlowingP  = globalFlowing[:, 1:kP]
        globalFlowingL  = globalFlowing[:, kP+1:kP+kL]

        # volume fractions
        phiMatrix   = vec(sum(globalMatrix, dims=2)) ./ densityP
        phiEnclosed = phiW .+ vec(sum(globalEnclosedP, dims=2)) ./ densityP .+
                      vec(sum(globalEnclosedL, dims=2)) ./ densityL
        phiBiofilm  = phiMatrix .+ phiEnclosed
        phiFlowing  = 1 .- phiBiofilm

        # clogging guard
        clog = findfirst(>(clogging_fraction), phiBiofilm)
        if clog !== nothing
            results.flag = "CLOGGED"
            results.simulation_data[:error] = (description="biofilm volume fraction exceeded clogging value",
                                               cell=clog, time=t)
            break
        end

        # local concentrations [X | S | quotients]
        localBiofilm = _local(globalMatrix ./ (phiBiofilm .+ _REALMIN),
                              globalEnclosedL ./ (phiBiofilm .+ _REALMIN),
                              kernel.num_idx, kernel.den_idx)
        localEnclosedX = globalEnclosedP ./ (phiEnclosed .+ _REALMIN)
        localEnclosed = _local(localEnclosedX,
                               globalEnclosedL ./ (phiEnclosed .+ _REALMIN),
                               kernel.num_idx, kernel.den_idx)
        localFlowingX = globalFlowingP ./ phiFlowing
        localFlowingS = globalFlowingL ./ phiFlowing
        localFlowing = _local(localFlowingX, localFlowingS, kernel.num_idx, kernel.den_idx)
        localFlowingXS = hcat(localFlowingX, localFlowingS)

        # attachment likelihoods
        attEnclosedFactor = phiMatrix ./ phiBiofilm
        replace!(attEnclosedFactor, NaN => 0.0)
        attFlowingFactor = (1 .- porosity_centers) .+ porosity_centers .* phiBiofilm

        # --- reaction terms: light ---
        globalParticles = globalMatrix .+ globalEnclosedP .+ globalFlowingP
        etaParticles = cumsum(vec(sum(attenuation_particles .* globalParticles, dims=2))) .* dz
        eta = eta_water .+ eta_sand .+ etaParticles
        light = f.light_irradiation(t)
        lightAttenuated = light .* exp.(-eta) ./ light_optimal
        lightEffective = lightAttenuated .* exp.(1 .- lightAttenuated)
        if any(light_dep)
            lightFactor[:, light_dep] = _light_factor_floor(lightEffective, min_light)
        end
        # Dark-switch reactions (Wolf2007 r6): K/(K + I_local), I_local in
        # optimal-intensity units (same normalization as lightAttenuated).
        for j in findall(inh_dep)
            lightFactor[:, j] .= light_inh[j] ./ (light_inh[j] .+ lightAttenuated)
        end
        # Complement reactions: 1 − Steele(I) (Steele ≤ 1, so the factor stays
        # in [0, 1]) — on in darkness, zero at optimal light.
        for j in findall(comp_dep)
            lightFactor[:, j] .= 1 .- lightEffective
        end

        # --- ecological + exchange reactions ---
        # Phase efficiencies scale each reaction per region (defaults 1.0). The
        # reaction rates are written into reused buffers and scaled in place.
        ecoBiofilm = _evaluate_reactions!(ecoBiofilmBuf, localBiofilm, kernel, phiBiofilm, muRates, lightFactor)
        ecoBiofilm .*= eff_biofilm
        ecoEnclosed = _evaluate_reactions!(ecoEnclosedBuf, localEnclosed, kernel, phiEnclosed, muRates, lightFactor)
        ecoEnclosed .*= eff_biofilm
        ecoFlowing = _evaluate_reactions!(ecoFlowingBuf, localFlowing, kernel, phiFlowing, muRates, lightFactor)
        ecoFlowing .*= eff_flowing

        ecoRxM  = ecoBiofilm * sigmaP'
        ecoRxPe = ecoEnclosed * sigmaP'
        ecoRxLe = (ecoBiofilm .+ ecoEnclosed) * sigmaL'
        ecoRxPf = ecoFlowing * sigmaP'
        ecoRxLf = ecoFlowing * sigmaL'

        attE = attEnclosedFactor .* globalEnclosedP .* attachment_rates
        # Flowing attachment splits into a sand term (scaled per particle by
        # `sand_attachment_factor`) and a biofilm term. With all factors 1.0 this
        # equals `attFlowingFactor .* globalFlowingP .* attachment_rates`.
        attFactorP = (1 .- porosity_centers) .* sand_factors .+ porosity_centers .* phiBiofilm
        attF = attFactorP .* globalFlowingP .* attachment_rates

        velFlowingCenters = 0.5 .* (velFlowing[2:end] .+ velFlowing[1:end-1])
        detM = model.detachment(velFlowingCenters) .* globalMatrix

        transP = (phiEnclosed ./ beta) .* (localFlowingX .- localEnclosedX) .* transport_particle_rates
        transL = (phiEnclosed ./ beta) .* (localFlowingS .- (globalEnclosedL ./ (phiEnclosed .+ _REALMIN))) .* transport_liquid_rates

        rxMatrix   = ecoRxM .+ attE .+ attF .- detM
        rxEnclP    = ecoRxPe .- attE .+ transP
        rxEnclL    = ecoRxLe .+ transL
        rxFlowP    = ecoRxPf .- attF .+ detM .- transP
        rxFlowL    = ecoRxLf .- transL

        rhsBiofilm = hcat(rxMatrix, rxEnclP, rxEnclL)
        rhsFlowing = hcat(rxFlowP, rxFlowL)
        rhsEnclWater = (beta .* phiBiofilm .- phiEnclosed) ./ tau
        # Implicit-osmosis option: the relaxation is (β·φB−φe)/τ = Aosm − kosm·φW
        # with kosm = (1−β)/τ, a stiff linear decay in φW. Integrating that term
        # with backward Euler damps the rate to rate/(1+dt·kosm) — unconditionally
        # stable, exact as dt→0, and zero at equilibrium (τ, β unchanged). This
        # lets us drop 1/τ from the CFL below so dt is no longer osmosis-bound.
        # Solver A runs on the previous dt (Diehl2025 §3.5); the φW update below
        # re-damps with the new dt.
        kosm = (1 - beta) / tau
        rhsEnclWaterA = implicit_osmosis ? rhsEnclWater ./ (1 .+ dt .* kosm) : rhsEnclWater

        # --- SOLVER A: biofilm velocity ---
        rhsBioVol = vec(sum(rxMatrix, dims=2)) ./ densityP .+
                    vec(sum(rxEnclP, dims=2)) ./ densityP .+
                    vec(sum(rxEnclL, dims=2)) ./ densityL .+ rhsEnclWaterA
        u = phiBiofilm[1:n0]
        uB = 0.5 .* (u[2:end] .+ u[1:end-1])               # length n0-1
        lambda = zeta_0 .* mobility(uB)
        Dmat = sparse(Dmob.I, Dmob.J, Dmob.V .* repeat(lambda, inner=4), 2n0, 2n0)
        lhsCH = CH0 - dt * (Smat + Dmat)
        rhsCH = vcat(u .+ dt .* rhsBioVol[1:n0], dpsi_fun(u))
        xCH = lhsCH \ rhsCH
        muCH = xCH[n0+1:end]                               # length n0
        velBiofilm[1:n0-1] = volume_avg_velocity[2:n0] .- zeta_0 .* (1 .- uB) .* diff(muCH) ./ dz

        # --- SOLVER B: transport update ---
        phiBdy = 0.5 .* (phiBiofilm[2:end] .+ phiBiofilm[1:end-1])
        velFlowing = vcat(volume_avg_velocity[1],
                          (volume_avg_velocity[2:end-1] .- velBiofilm .* phiBdy) ./ (1 .- phiBdy),
                          volume_avg_velocity[end])

        # --- time adaptivity (CFL) ---
        # Per-step CFL bound; dt grows by at most (1 + tol) toward
        # cfl_factor/max(w). Weights per region (matrix, enclosed P, enclosed L,
        # flowing P, flowing L): w_v advection, w_a dispersion, w_b exchange,
        # w_s ecological source. Ported from the MATLAB adaptive path.
        if is_adaptive
            vbmax = (1 + adaptive_velocity_factor) * maximum(abs.(velBiofilm))
            vfmax = (1 + adaptive_velocity_factor) * maximum(abs.(velFlowing))
            maxPhib = maximum(phiBiofilm)
            phie_f_max = maximum(phiEnclosed ./ phiFlowing)
            det_vf = model.detachment(velFlowingCenters)

            # region-local particle (X) and liquid (S) concentrations
            Xb = globalMatrix ./ (phiBiofilm .+ _REALMIN)
            Sb = globalEnclosedL ./ (phiBiofilm .+ _REALMIN)
            Xe = localEnclosedX
            Se = globalEnclosedL ./ (phiEnclosed .+ _REALMIN)
            Xf = localFlowingX
            Sf = localFlowingS

            # liquid-consumption bound (growth reactions) per region
            L_b = _cfl_liquid_bound(Xb, Sb, sigmaL, muRates, K_CFL, kL)
            L_e = _cfl_liquid_bound(Xe, Se, sigmaL, muRates, K_CFL, kL)
            L_f = _cfl_liquid_bound(Xf, Sf, sigmaL, muRates, K_CFL, kL)

            # hydrolysis quotient monod (POM/HET = particle columns 3, 1)
            maxXi(X) = _cfl_hydrolysis_bound(X, K_Hyd)
            s35mu5 = abs(sigmaP[3, 5]) * muRates[5]
            ws_t0 = max(abs(sigmaP[1, 3]) * muRates[3], abs(sigmaP[2, 4]) * muRates[4])

            w_v = [2vbmax, 2vbmax, 2vbmax, 2vfmax, 2vfmax]
            w_a = [0.0, 0.0, 0.0,
                   2vfmax * alphaP * (1 + 1 / (1 - maxPhib)),
                   2vfmax * alphaL * (1 + 1 / (1 - maxPhib))]
            w_b = [maximum(det_vf),
                   maximum(attachment_rates) * maximum(attEnclosedFactor) + maximum(transport_particle_rates) / beta,
                   implicit_osmosis ? maximum(transport_liquid_rates) / beta :
                       max(maximum(transport_liquid_rates) / beta, 1 / tau),
                   maximum(attachment_rates) * maximum(attFlowingFactor) + maximum(transport_particle_rates) / beta * phie_f_max,
                   maximum(transport_liquid_rates) / beta * phie_f_max]
            w_s = [max(ws_t0, s35mu5 * maxXi(Xb)),
                   max(ws_t0, s35mu5 * maxXi(Xe)),
                   L_b + L_e,
                   max(ws_t0, s35mu5 * maxXi(Xf)),
                   L_f]

            dt_CFL = cfl_factor / maximum(w_v ./ dz .+ w_a ./ dz^2 .+ w_b .+ w_s)
            dt = min(dt_CFL, (1 + adaptive_time_tolerance) * dt, adaptive_max_dt)
        end

        zrowB = zeros(1, 2kP + kL)
        zrowF = zeros(1, kP + kL)

        # biofilm flux (upwind on velBiofilm)
        fBio = vcat(zrowB,
                    globalBiofilm[1:end-1, :] .* max.(0.0, velBiofilm) .+
                    globalBiofilm[2:end, :] .* min.(0.0, velBiofilm),
                    zrowB)
        fBioIn  = porosity_boundaries[1:end-1] .* fBio[1:end-1, :]
        fBioOut = porosity_boundaries[2:end] .* fBio[2:end, :]

        # flowing: convection (upwind on velFlowing) minus dispersion
        localGrad = diff(localFlowingXS[1:end-1, :]; dims=1)           # (N-2)×(kP+kL)
        dispStrength = abs.(velFlowing[2:end-2]) .* (1 .- phiBdy[1:end-1])
        dispFlux = vcat(zrowF,
                        (dispStrength .* localGrad .* alpha) ./ dz,
                        zrowF, zrowF)
        convFlux = vcat(volume_avg_velocity[1] .* globalConcInflow,
                        globalFlowing[1:end-1, :] .* max.(0.0, velFlowing[2:end-1]) .+
                        globalFlowing[2:end, :] .* min.(0.0, velFlowing[2:end-1]),
                        volume_avg_velocity[end] .* globalFlowing[end:end, :])
        fFlow = convFlux .- dispFlux
        fFlowIn  = porosity_boundaries[1:end-1] .* fFlow[1:end-1, :]
        fFlowOut = porosity_boundaries[2:end] .* fFlow[2:end, :]

        # enclosed water flux (upwind on velBiofilm)
        fWat = vcat(0.0,
                    phiW[1:end-1] .* max.(0.0, velBiofilm) .+ phiW[2:end] .* min.(0.0, velBiofilm),
                    0.0)
        fWatIn  = porosity_boundaries[1:end-1] .* fWat[1:end-1]
        fWatOut = porosity_boundaries[2:end] .* fWat[2:end]

        # update cell values
        globalBiofilm = globalBiofilm .+ (dt / dz) .* (fBioIn .- fBioOut) ./ porosity_centers .+ dt .* rhsBiofilm
        globalFlowing = globalFlowing .+ (dt / dz) .* (fFlowIn .- fFlowOut) ./ porosity_centers .+ dt .* rhsFlowing
        rhsEnclWaterB = implicit_osmosis ? rhsEnclWater ./ (1 .+ dt .* kosm) : rhsEnclWater
        phiW = phiW .+ (dt / dz) .* (fWatIn .- fWatOut) ./ porosity_centers .+ dt .* rhsEnclWaterB

        # Underflow floor: a sink acting on an exactly-zero pool (r6 consuming
        # PG before any is stored; POM hydrolysis) leaves O(realmin) negatives
        # from the Monod regularizer that would trip the strict guard below.
        # Clamp only denormal-scale negatives; genuine instabilities overshoot
        # far beyond -1e-20.
        @. globalBiofilm = ifelse(-1e-20 < globalBiofilm < 0, 0.0, globalBiofilm)
        @. globalFlowing = ifelse(-1e-20 < globalFlowing < 0, 0.0, globalFlowing)

        # negativity guards
        if any(x -> x < 0 || isnan(x), globalBiofilm)
            results.flag = "BIOFILM"
            results.simulation_data[:error] = (description="unphysical concentration in biofilm", time=t)
            break
        elseif any(x -> x < 0 || isnan(x), globalFlowing)
            results.flag = "FLOWING"
            results.simulation_data[:error] = (description="unphysical concentration in flowing suspension", time=t)
            break
        end

        # advance + frame capture
        t += dt
        push!(step_times, t)
        if counter <= numFrames && t >= timeSnap[counter]
            timeFrames[counter] = t
            concFramesBiofilm[:, counter, :] = globalBiofilm
            concFramesFlowing[:, counter, :] = globalFlowing
            concFramesWater[:, counter] = densityL .* phiW
            velFramesBiofilm[:, counter] = velBiofilm
            velFramesFlowing[:, counter] = velFlowing[2:end-1]
            concFramesReaction[:, counter, :] = ecoBiofilm
            counter += 1
            quiet || @info "simulate" t
        end
    end

    # ---- VII. outputs ------------------------------------------------------
    results.simulation_data[:step_times] = step_times
    results.frames[:time] = timeFrames
    results.frames[:concentration_biofilm] = concFramesBiofilm
    results.frames[:concentration_flowing] = concFramesFlowing
    results.frames[:concentration_water] = concFramesWater
    results.frames[:velocity_biofilm] = velFramesBiofilm
    results.frames[:velocity_flowing] = velFramesFlowing
    results.frames[:reaction_rates] = concFramesReaction
    results.time_final = t
    results.simulation_data[:time_final] = t
    # For adaptive runs keep the "adaptive" marker (record the last dt separately);
    # for fixed runs store the constant dt.
    if is_adaptive
        results.simulation_data[:time_step] = "adaptive"
        results.simulation_data[:final_dt] = dt
    else
        results.simulation_data[:time_step] = dt
    end
    return results
end
