# Render the Log-OAT sensitivity outputs (per Logarithmic_OAT_Sensitivity_SSF_model.pdf):
#   * a horizontal bar chart of the RMS log-sensitivity I_i (OK params),
#   * the L+_i − L0 / L−_i − L0 delta curves for the top parameters.
# Run with the plotenv (MPCSSF dev + CairoMakie):
#   julia --project=<plotenv> julia/analysis/render_log_oat.jl
using CairoMakie, DelimitedFiles

dir = joinpath(@__DIR__, "results", "log_oat")
M, hdr = readdlm(joinpath(dir, "measures.csv"), ','; header=true)
names = String.(M[:, 1]); blocks = String.(M[:, 2])
Irms = Float64.(M[:, 3]); clog = Bool.(M[:, 10] .== true)

# ---- I_i ranking bar chart (OK params, descending) --------------------------
ok = .!clog
order = sortperm(Irms[ok]; rev=false)          # ascending for horizontal bars
nm = names[ok][order]; val = Irms[ok][order]
fig = Figure(size=(720, 640))
ax = Axis(fig[1, 1]; xscale=log10, xlabel="RMS log-sensitivity  Iᵢ  (log scale)",
          title="Log-OAT sensitivity of pathogen log-removal (mature + pulse)",
          yticks=(1:length(nm), nm))
barplot!(ax, 1:length(nm), max.(val, 1e-7); direction=:x, color=:steelblue)
save(joinpath(dir, "ranking_Irms.png"), fig)

# ---- delta curves for the top-6 OK params + note clog-driving ---------------
top = names[ok][sortperm(Irms[ok]; rev=true)][1:min(6, count(ok))]
fig2 = Figure(size=(900, 620))
for (i, p) in enumerate(top)
    C, _ = readdlm(joinpath(dir, "curves_$(p).csv"), ','; header=true)
    t = Float64.(C[:, 1]); dLp = Float64.(C[:, 2]); dLm = Float64.(C[:, 3])
    r, c = fldmod1(i, 3)
    ax = Axis(fig2[r, c]; title=p, xlabel="t (day)", ylabel="Lᵢ − L₀")
    lines!(ax, t, dLp; color=:firebrick, label="×2")
    lines!(ax, t, dLm; color=:seagreen, label="×½")
    hlines!(ax, [0.0]; color=(:gray, 0.5), linestyle=:dash)
    i == 1 && axislegend(ax; position=:rb)
end
Label(fig2[0, :], "Response curves Lᵢ(t) − L₀(t) for the top-6 parameters"; fontsize=16)
save(joinpath(dir, "delta_curves_top6.png"), fig2)
println("wrote ranking_Irms.png and delta_curves_top6.png to ", dir)
