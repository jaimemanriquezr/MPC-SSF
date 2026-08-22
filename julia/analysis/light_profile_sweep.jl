# Light–biofilm-profile sweep: vary irradiance amplitude at fixed geometry and
# report where biofilm accumulates.
#
# Two QoIs, per Jaime: (1) the depth of the peak biofilm volume fraction, and
# (2) biofilm mass above z = 0 (the supernatant + roughness layer).
#
# Geometry is fixed at sand_roughness = 20 mm so the near-optimal light band sits
# ~8 mm above the nominal sand surface, in water that is ~36% solid and therefore
# supports bare-sand attachment. See
# .claude/plans/2026-07-29-light-biofilm-profile.md.
#
# FIXED time stepping, not adaptive: the adaptive CFL block hardcodes the Lund
# structure (it indexes sigmaParticles[3,5] and muRates[3:5]), so it throws a
# BoundsError on this 2-particle/4-reaction model. Costs nothing — measured legacy
# scaling is exactly linear in cell count, which shows dt is pinned by the cap
# rather than by a CFL bound.
#
# Run: julia --project=. analysis/light_profile_sweep.jl [amplitude] [tsim] [ncells]

using SSF, Printf, DelimitedFiles, Statistics
include(joinpath(@__DIR__, "light_profile_model.jl"))

const DELTA = 20e-3   # sand_roughness, held fixed across the sweep

"""
    profile_qois(r, f) -> NamedTuple

`z_peak` is the depth (m) of the maximum biofilm volume fraction in the final
frame; `mass_above` is ∫ porosity·φ_b dz over z < 0; `mass_total` the same over the
whole domain. Frames are truncated at `time_final` because an early break leaves
preallocated tail frames as zeros at time 0.0.
"""
function profile_qois(r, f)
    vf = get_volume_fractions(r).biofilm
    z  = depths(r); dz = SSF.gridsize(f)
    poro = SSF.computeporosity(f, z)
    ts = times(r)
    # Frames are preallocated and filled left-to-right; on an early break the tail
    # stays zeros AT TIME 0.0. Frame 1 is legitimately t=0, so count from frame 2.
    nvalid = 1 + count(>(0), @view ts[2:end])
    keep = 1:nvalid
    last_frame = nvalid
    phib = vf[:, last_frame]
    above = z .< 0
    return (z_peak = z[argmax(phib)],
            phib_peak = maximum(phib),
            mass_above = sum(poro[above] .* phib[above]) * dz,
            mass_total = sum(poro .* phib) * dz,
            phib_max_above = maximum(phib[above]),
            nframes_valid = length(keep),
            profile = phib, z = z, times = ts[keep],
            mass_above_traj = [sum(poro[above] .* vf[above, k]) * dz for k in keep],
            mass_total_traj = [sum(poro .* vf[:, k]) * dz for k in keep])
end

function run_one(amplitude; tsim=7.0, ncells=500, nframes=60, dt=3e-6, quiet=true)
    m = light_profile_model()
    f = light_filter(amplitude=amplitude, ncells=ncells, delta=DELTA)
    t0 = time()
    r = simulate(State(f, m); inflow_concentrations=INFLUENT_LIGHT, simulation_time=tsim,
                 time_step=dt, n_frames=nframes, quiet=quiet, implicit_osmosis=true)
    wall = time() - t0
    q = profile_qois(r, f)
    @printf("A=%-7.4g flag=%-8s t_final=%.4f wall=%.0fs (%.0f s/sim-day)  z_peak=%+8.2f mm  phib_peak=%.4g  M_above=%.4g  M_total=%.4g  max_phib(z<0)=%.4g\n",
            amplitude, string(r.flag), r.time_final, wall, wall/max(r.time_final,1e-9),
            q.z_peak*1000, q.phib_peak, q.mass_above, q.mass_total, q.phib_max_above)
    return (; amplitude, flag=string(r.flag), t_final=r.time_final, wall, q...)
end

function main(amplitude, tsim, ncells, dt)
    outdir = joinpath(@__DIR__, "results", "light_profile"); isdir(outdir) || mkpath(outdir)
    @info "light profile" amplitude tsim ncells dt delta_mm=DELTA*1e3
    res = run_one(amplitude; tsim=tsim, ncells=ncells, dt=dt, quiet=true)
    tag = @sprintf("A%.4g", amplitude)
    writedlm(joinpath(outdir, "profile_$tag.csv"),
             vcat(["z" "phi_b"], hcat(res.z, res.profile)), ',')
    writedlm(joinpath(outdir, "traj_$tag.csv"),
             vcat(["t" "mass_above" "mass_total"],
                  hcat(res.times, res.mass_above_traj, res.mass_total_traj)), ',')
    writedlm(joinpath(outdir, "qoi_$tag.csv"),
             ["amplitude" "flag" "t_final" "z_peak_m" "phib_peak" "mass_above" "mass_total" "max_phib_above" "wall_s";
              res.amplitude res.flag res.t_final res.z_peak res.phib_peak res.mass_above res.mass_total res.phib_max_above res.wall], ',')
    println("wrote ", outdir, " (", tag, ")")
end

if abspath(PROGRAM_FILE) == @__FILE__
    A      = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.8
    tsim   = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 7.0
    ncells = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 500
    # dt default 3e-7, NOT run_proxy's 3e-6. Measured: 3e-6 trips the biofilm
    # negativity guard at t=0.2323 at ncells=500 -- and so does plain modelLund
    # (t=0.2322), so this is the mesh/dt combination, not the reduced model.
    # 1e-6 survived 0.3 d; 3e-7 is the margin for a 7-day run as biofilm grows
    # and tightens the CFL.
    dt     = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 3e-7
    main(A, tsim, ncells, dt)
end
