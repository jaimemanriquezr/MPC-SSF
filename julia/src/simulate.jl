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

# Port of the evaluateReactions subfunction. `phi` is per-cell; returns nCells×nRx.
function _evaluate_reactions(local_, kernel, phi, muRates, lightFactor)
    nCells = size(local_, 1)
    nRx = length(muRates)
    rx = zeros(nCells, nRx)
    for i in 1:nRx
        cols = kernel.monod_cols[i]
        if isempty(cols)
            monod = ones(nCells)
        else
            sub = local_[:, cols]
            Ki = reshape(kernel.monod_K[i], 1, :)
            terms = (sub .+ _REALMIN) ./ (Ki .+ sub .+ _REALMIN)
            monod = vec(minimum(terms, dims=2))
        end
        prod_ = @view local_[:, kernel.order_col[i]]
        @views @. rx[:, i] = phi * muRates[i] * lightFactor[:, i] * monod * prod_
    end
    return rx
end

_meanval(x) = sum(x) / length(x)

"""
    simulate(state::State; inflow_concentrations=nothing, simulation_time=1.0,
             time_step=1e-5, n_frames=200, clogging_fraction=0.99,
             upwinded=false, quiet=false) -> Results

Run the slow-sand-filtration simulation from `state` with a fixed time step
(parity with the current MATLAB `simulate`). `inflow_concentrations` may be
`nothing` (zero inflow), a vector of length `kP+kL`, or a function `t -> vector`.
Returns a [`Results`](@ref) whose `flag` is `"OK"`, `"CLOGGED"`, `"BIOFILM"`, or
`"FLOWING"`.
"""
function simulate(state::State;
                  inflow_concentrations=nothing,
                  simulation_time::Real=1.0,
                  time_step::Real=1e-5,
                  n_frames::Integer=200,
                  clogging_fraction::Real=0.99,
                  upwinded::Bool=false,
                  quiet::Bool=false)
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
    light_optimal = isempty(model.reactions) ? 1.0 :
                    maximum(r.optimal_light_factor for r in model.reactions)
    attenuation_particles = reshape([p.attenuation for p in P], 1, kP)
    light_dep = [r.is_light_dependent for r in model.reactions]
    min_light = [r.minimum_light_factor for r in model.reactions][light_dep]
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
    dt = float(time_step)
    t = t0
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
        lightFactor = ones(N, nRx)
        if any(light_dep)
            tmp = reshape(min_light, 1, :) .+ lightEffective       # N × nLightDep
            lightFactor[:, light_dep] = max.(0.0, tmp)
        end

        # --- ecological + exchange reactions ---
        ecoBiofilm = _evaluate_reactions(localBiofilm, kernel, phiBiofilm, muRates, lightFactor)
        ecoEnclosed = _evaluate_reactions(localEnclosed, kernel, phiEnclosed, muRates, lightFactor)
        ecoFlowing = _evaluate_reactions(localFlowing, kernel, phiFlowing, muRates, lightFactor)

        ecoRxM  = ecoBiofilm * sigmaP'
        ecoRxPe = ecoEnclosed * sigmaP'
        ecoRxLe = (ecoBiofilm .+ ecoEnclosed) * sigmaL'
        ecoRxPf = ecoFlowing * sigmaP'
        ecoRxLf = ecoFlowing * sigmaL'

        attE = attEnclosedFactor .* globalEnclosedP .* attachment_rates
        attF = attFlowingFactor .* globalFlowingP .* attachment_rates

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

        # --- SOLVER A: biofilm velocity ---
        rhsBioVol = vec(sum(rxMatrix, dims=2)) ./ densityP .+
                    vec(sum(rxEnclP, dims=2)) ./ densityP .+
                    vec(sum(rxEnclL, dims=2)) ./ densityL .+ rhsEnclWater
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
        phiW = phiW .+ (dt / dz) .* (fWatIn .- fWatOut) ./ porosity_centers .+ dt .* rhsEnclWater

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
        if counter <= numFrames && t >= timeSnap[counter]
            timeFrames[counter] = t
            concFramesBiofilm[:, counter, :] = globalBiofilm
            concFramesFlowing[:, counter, :] = globalFlowing
            concFramesWater[:, counter] = densityL .* phiW
            velFramesBiofilm[:, counter] = velBiofilm
            velFramesFlowing[:, counter] = velFlowing[2:end-1]
            counter += 1
            quiet || @info "simulate" t
        end
    end

    # ---- VII. outputs ------------------------------------------------------
    results.frames[:time] = timeFrames
    results.frames[:concentration_biofilm] = concFramesBiofilm
    results.frames[:concentration_flowing] = concFramesFlowing
    results.frames[:concentration_water] = concFramesWater
    results.frames[:velocity_biofilm] = velFramesBiofilm
    results.frames[:velocity_flowing] = velFramesFlowing
    results.simulation_data[:time_final] = t
    results.simulation_data[:time_step] = dt
    return results
end
