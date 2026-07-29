# Variance-based (Sobol) global sensitivity — UQ plan stage 3. Retains the most
# influential NON-clogging parameters from the Log-OAT screen and estimates
# first-order Sᵢ and total-effect S_Tᵢ indices for the pathogen-removal QoI
# (mean log-removal L̄ over the post-disturbance window, mature+pulse scenario).
#
# Design: Saltelli sampling of two independent N×k unit matrices A, B and the k
# hybrid matrices AB⁽ⁱ⁾ (A with column i from B) → N·(k+2) model runs. Estimators:
#   Sᵢ   = mean(y_B · (y_ABi − y_A)) / Var(Y)             (Saltelli 2010)
#   S_Tᵢ = mean((y_A − y_ABi)²) / (2·Var(Y))              (Jansen 1999)
# Bootstrap CIs by resampling the N realizations (no extra model runs).
#
# The clog-driving params (velocity, beta_porosity, zeta_0) are EXCLUDED — their
# ranges cause filter failure (undefined QoI); they were reported separately by
# the OAT. Run on the fast proxy so N·(k+2) is workstation-scale.
#
# Run: julia --project=julia julia/analysis/sobol_sensitivity.jl [N] [ncells] [tpost]
# Writes results/sobol/sobol_indices.csv and prints ranked Sᵢ / S_Tᵢ with CIs.

using MPCSSF, Printf, DelimitedFiles, Statistics, Random
include(joinpath(@__DIR__, "log_oat_sensitivity.jl"))   # builders, run_proxy, build_mature, make_inflow, TD, PAT_IN

struct SP; name::String; lo::Float64; hi::Float64; scale::Symbol; apply::Function; end
_val(p::SP, u) = p.scale === :log ? p.lo * (p.hi / p.lo)^u : p.lo + (p.hi - p.lo) * u

# retained parameters + ranges (cited from Campos2006 / Schijven2013, as in the OAT)
const SPARAMS = SP[
    SP("attach_sand",  180.0,  1640.0, :log, set_attach_sand()),                  # Lund b_sand ±3×
    SP("sand_pathogen", 5e-3,   0.71,  :log, set_sand_pat()),                     # Schijven α
    SP("dispersivity",  3.6e-3, 3.6e-2,:log, set_disp()),                         # Schijven T1
    SP("temperature",   3.0,    25.0,  :lin, set_temp()),                         # seasonal
    SP("bacterivory",   2.0,    20.0,  :log, set_rate(["Bacterivory"])),          # Manriquez B.4
    SP("theta_death",   0.03,   0.12,  :lin, set_theta(["Heterotroph death","Phototroph death"])),  # θ−1
    SP("transport_P",   1.64,   16.4,  :log, set_transport_P()),                  # Lund
    SP("influent_PAT",  0.3,    3.0,   :log, (f,m,v)->(f,m)),                      # multiplier (special)
]

# QoI: mean log-removal over the post-disturbance window (mature+pulse)
function eval_qoi(u::AbstractVector, ms, m0, ncells, tpost, nframes)
    f2 = addgridpoints(SandFilter(temperature=15.0), ncells); m2 = m0; pat_mult = 1.0
    for (p, ui) in zip(SPARAMS, u)
        v = _val(p, ui)
        p.name == "influent_PAT" ? (pat_mult = v) : ((f2, m2) = p.apply(f2, m2, v))
    end
    s = State(f2, m2, 0.0, deepcopy(ms.global_concentration),
              copy(ms.enclosed_water_volume), deepcopy(ms.velocity))
    infl, Cref = make_inflow(:pulse, pat_mult)
    r = run_proxy(s; simulation_time=tpost, inflow_concentrations=infl, n_frames=nframes, quiet=true)
    ts = times(r); cout = max.(concentration(r, "PAT", :flowing)[end, :], 1e-30)
    L = log10.(Cref ./ cout); win = (ts .>= TD) .& (ts .<= r.time_final + 1e-9)
    return (mean(L[win]), string(r.flag))
end

# bootstrap CI for an estimator functional over the N realizations
function _bootstrap(f, N, rng; B=500)
    vals = [f(rand(rng, 1:N, N)) for _ in 1:B]
    return quantile(vals, 0.05), quantile(vals, 0.95)
end

function main(; N=48, ncells=25, tpost=1.0, nframes=30, tmature=3.0, seed=20260722)
    rng = MersenneTwister(seed)
    k = length(SPARAMS)
    outdir = joinpath(@__DIR__, "results", "sobol"); isdir(outdir) || mkpath(outdir)
    @info "Sobol" N k runs=N*(k+2) ncells tpost

    m0 = pathogen_model(); ms = build_mature(m0, ncells, tmature)
    A = rand(rng, N, k); B = rand(rng, N, k)
    ev(row) = eval_qoi(row, ms, m0, ncells, tpost, nframes)

    # Each row is an independent simulation (own filter, own State, only reads
    # `ms`), so the design matrix parallelises across threads. Needed to make a
    # useful N reachable: N*(k+2) runs at ~1.7 min each is days on one core.
    # nclog is atomic because threads increment it concurrently.
    nclog = Threads.Atomic{Int}(0)
    function run_matrix(M)
        y = Vector{Float64}(undef, size(M, 1))
        Threads.@threads for i in 1:size(M, 1)
            yi, fl = ev(@view M[i, :])
            fl == "OK" || Threads.atomic_add!(nclog, 1)
            y[i] = yi
        end
        return y
    end
    yA = run_matrix(A); @printf("A done (%d clog)\n", nclog[])
    yB = run_matrix(B); @printf("B done (%d clog cum)\n", nclog[])

    varY = var(vcat(yA, yB))
    rows = Vector{Any}[["param","S_i","S_i_lo","S_i_hi","S_Ti","S_Ti_lo","S_Ti_hi"]]
    results = Tuple{String,Float64,Float64,Float64,Float64,Float64}[]
    for i in 1:k
        ABi = copy(A); ABi[:, i] = B[:, i]
        yABi = run_matrix(ABi)
        Si  = mean(yB .* (yABi .- yA)) / varY
        STi = mean((yA .- yABi) .^ 2) / (2 * varY)
        # bootstrap CIs (resample realizations)
        rb = MersenneTwister(seed + i)
        silo, sihi = _bootstrap(idx -> mean(yB[idx] .* (yABi[idx] .- yA[idx])) / var(vcat(yA[idx], yB[idx])), N, rb)
        stlo, sthi = _bootstrap(idx -> mean((yA[idx] .- yABi[idx]) .^ 2) / (2 * var(vcat(yA[idx], yB[idx]))), N, rb)
        push!(rows, [SPARAMS[i].name, Si, silo, sihi, STi, stlo, sthi])
        push!(results, (SPARAMS[i].name, Si, STi, silo, sihi, stlo))
        @printf("  %-14s  S_i=%+.3f [%.3f,%.3f]   S_Ti=%.3f [%.3f,%.3f]\n",
                SPARAMS[i].name, Si, silo, sihi, STi, stlo, sthi)
    end
    writedlm(joinpath(outdir, "sobol_indices.csv"), rows, ',')
    @printf("\nVar(Y)=%.4f  total runs=%d  (%d non-OK)\n", varY, N*(k+2), nclog[])
    println("\n── Ranked by total-effect S_Ti ──")
    for (nm, si, sti) in sort([(r[1], r[2], r[3]) for r in results]; by=x->x[3], rev=true)
        @printf("  S_Ti=%.3f  S_i=%+.3f  interaction=%.3f  %s\n", sti, si, sti - si, nm)
    end
    println("\nwrote ", joinpath(outdir, "sobol_indices.csv"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    N      = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 48
    ncells = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 25
    tpost  = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1.0
    main(; N=N, ncells=ncells, tpost=tpost)
end
