# Numerical-uncertainty / mesh-convergence study (UQ plan step 1): repeat the
# mature+pulse baseline at several grid resolutions and check that the sensitivity
# QoIs (mean/min log-removal, biofilm mass) converge. The requirement is that the
# numerical change in each QoI between successive meshes is small relative to the
# parameter-induced spread found by the OAT (top I_rms ~0.5 in L), i.e. ≲1–2%.
#
# Run: julia --project=julia julia/analysis/mesh_convergence.jl [tmature] [tpost]
# Writes results/mesh_convergence.csv and prints a table with successive rel-changes.

using SSF, Printf, DelimitedFiles, Statistics
include(joinpath(@__DIR__, "log_oat_sensitivity.jl"))   # pathogen_model, run_proxy, challenge_inflow, INFLUENT_BASE, TD

function qoi_at_mesh(ncells, tmature, tpost, nframes)
    m0 = pathogen_model()
    f = addgridpoints(SandFilter(temperature=15.0), ncells)
    rmat = run_proxy(State(f, m0); simulation_time=tmature, inflow_concentrations=INFLUENT_BASE,
                     n_frames=5, quiet=true)
    z = depths(rmat); dz = gridsize(f); poro = computeporosity(f, z)
    phib = get_volume_fractions(rmat).biofilm[:, end]
    Mb = sum(poro .* phib .* dz)                              # biofilm mass ∫φ_b dz
    zb = sum(z .* phib) / (sum(phib) + 1e-30)                 # centroid
    ms = final_state(rmat)
    s = State(addgridpoints(SandFilter(temperature=15.0), ncells), m0, 0.0,
              deepcopy(ms.global_concentration), copy(ms.enclosed_water_volume), deepcopy(ms.velocity))
    infl, Cref = challenge_inflow(1.0)
    r = run_proxy(s; simulation_time=tpost, inflow_concentrations=infl, n_frames=nframes, quiet=true)
    ts = times(r); cout = max.(concentration(r, "PAT", :flowing)[end, :], 1e-30); L = log10.(Cref ./ cout)
    win = (ts .>= TD) .& (ts .<= r.time_final + 1e-9)
    return (; Lmean=mean(L[win]), Lmin=minimum(L[win]), Mb, zb,
            matflag=string(rmat.flag), flag=string(r.flag), phibmax=maximum(phib))
end

function main(; meshes=[15, 20, 30, 40, 50], tmature=3.0, tpost=1.5, nframes=60)
    outdir = joinpath(@__DIR__, "results"); isdir(outdir) || mkpath(outdir)
    @info "mesh-convergence" meshes tmature tpost
    Q = Any[]
    for nc in meshes
        q = qoi_at_mesh(nc, tmature, tpost, nframes)
        push!(Q, (nc, q))
        @printf("  ncells=%3d  Lmean=%.4f  Lmin=%.4f  Mb=%.4e  zb=%.4f  phibmax=%.3f  (%s/%s)\n",
                nc, q.Lmean, q.Lmin, q.Mb, q.zb, q.phibmax, q.matflag, q.flag)
    end
    # successive relative change (vs finest mesh as reference)
    ref = Q[end][2]
    rows = Vector{Any}[["ncells", "Lmean", "Lmin", "Mb", "zb", "relerr_Lmean", "relerr_Lmin", "relerr_Mb", "flag"]]
    println("\nRelative error vs finest mesh (ncells=$(meshes[end])):")
    for (nc, q) in Q
        rL = abs(q.Lmean - ref.Lmean) / (abs(ref.Lmean) + 1e-30)
        rm = abs(q.Lmin - ref.Lmin) / (abs(ref.Lmin) + 1e-30)
        rM = abs(q.Mb - ref.Mb) / (abs(ref.Mb) + 1e-30)
        push!(rows, [nc, q.Lmean, q.Lmin, q.Mb, q.zb, rL, rm, rM, q.flag])
        @printf("  ncells=%3d  ΔLmean=%.2f%%  ΔLmin=%.2f%%  ΔMb=%.2f%%\n", nc, 100rL, 100rm, 100rM)
    end
    writedlm(joinpath(outdir, "mesh_convergence.csv"), rows, ',')
    println("\nGuidance: numerical Δ should be ≲1–2% of the parameter-induced spread")
    println("(top OAT I_rms≈0.5 → ×2 changes L by ~0.37; so aim ΔLmean ≲ ~0.005 abs).")
    println("wrote ", joinpath(outdir, "mesh_convergence.csv"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    tmat = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 3.0
    tpost = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.5
    main(; tmature=tmat, tpost=tpost)
end
