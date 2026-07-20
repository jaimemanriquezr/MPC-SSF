# One-at-a-time (OAT) local sensitivity analysis for the pathogen slow-sand-filter
# model, for the Manriquez2026 revision (reviewers asked which parameters most
# influence the results). Scaffolding: a base case, a set of parameter
# perturbations grouped by the quantities the reviewers named
# (attachment/detachment, biofilm excess velocity, light attenuation, reaction
# rates, PAT inactivation/bacterivory), headline output metrics, and a driver
# that runs the base + each perturbation and reports relative sensitivities.
#
# Run:  julia --project=julia julia/analysis/sensitivity.jl
# Writes: julia/analysis/sensitivity_results.csv  and prints a ranked table.
#
# NOTE this is a *local* OAT sweep around one base case — it shows direction and
# rough magnitude of influence, not global/interaction effects. Extend `PERTURB`
# and `RUN` (horizon, seed, base case) as the study needs.

using MPCSSF
using Printf, DelimitedFiles

# ---- base run configuration -------------------------------------------------
const NCELL  = 40
const TSIM   = 0.3          # medium horizon (days): schmutzdecke operating, still transient
const DT     = 2e-5
const NFRAMES = 20
const INFLOW = Float64[1e-2, 1e-2, 0.0, 1e-3, 1e-2, 1e-2, 1e-3, 0.0, 1e-2]  # HET PHO POM PAT|O2 IC NH4 HPO4 DOM

# Seed a mature schmutzdecke near osmotic equilibrium (phi_e ≈ 99·phi_matrix), so
# an established biofilm operates from t=0 and pathogen removal is active. Same
# seed for every run, so a perturbation's *effect* is isolated even though the
# absolute baseline carries the seed.
function seed_mature!(s, z)
    for (i, zi) in enumerate(z)
        if zi >= 0
            d = exp(-zi / 0.15)
            s.global_concentration.matrix[i, 1] = 0.004 * d
            s.global_concentration.matrix[i, 2] = 0.002 * d
            s.global_concentration.matrix[i, 3] = 0.001 * d
            s.global_concentration.matrix[i, 4] = 2e-4 * d
            s.enclosed_water_volume[i] = 0.594 * d
            for l in 1:5
                s.global_concentration.enclosed_liquids[i, l] = 1e-3
            end
        end
    end
    return s
end

base_filter() = addgridpoints(SandFilter(temperature=20.0), NCELL)
base_model()  = modelPathogen()

# ---- immutable-struct "remake" helpers (perturb one field, keep the rest) ----
remake(rx::Reaction; kw...) = Reaction(; name=rx.name,
    (; nominal_rate=rx.nominal_rate, temperature_correction_factor=rx.temperature_correction_factor,
       order=rx.order, half_saturation_constants=rx.half_saturation_constants,
       stoichiometric_coefficients=rx.stoichiometric_coefficients,
       efficiency_biofilm=rx.efficiency_biofilm, efficiency_flowing=rx.efficiency_flowing,
       is_light_dependent=rx.is_light_dependent, minimum_light_factor=rx.minimum_light_factor,
       optimal_light_factor=rx.optimal_light_factor, kw...)...)

remake(p::Particle; kw...) = Particle(; name=p.name,
    (; density=p.density, dispersivity=p.dispersivity, transport_rate=p.transport_rate,
       attachment_matrix=p.attachment_matrix, attachment_sand=p.attachment_sand,
       attenuation=p.attenuation, sand_attachment_factor=p.sand_attachment_factor, kw...)...)

remake_model(m::Model; kw...) = Model(;
    (; components=m.components, reactions=m.reactions, cohesion_submodel=m.cohesion_submodel,
       water_density=m.water_density, biofilm_porosity=m.biofilm_porosity,
       osmosis_rate=m.osmosis_rate, detachment=m.detachment, kw...)...)

# model perturbation builders (m -> m')
scale_attach(k) = m -> remake_model(m; components=Component[
    c isa Particle ? remake(c; attachment_sand=c.attachment_sand*k, attachment_matrix=c.attachment_matrix*k) : c
    for c in m.components])
set_sand_pat(v) = m -> remake_model(m; components=Component[
    (c isa Particle && c.name=="PAT") ? remake(c; sand_attachment_factor=v) : c for c in m.components])
scale_detach(k) = m -> remake_model(m; detachment = v -> k .* m.detachment(v))
scale_zeta0(k)  = m -> remake_model(m; cohesion_submodel=CahnHilliardModel(
    kappa=m.cohesion_submodel.kappa, zeta_0=m.cohesion_submodel.zeta_0*k, zeta_1=m.cohesion_submodel.zeta_1))
scale_rx(idxs, k) = m -> remake_model(m; reactions=Reaction[
    (i in idxs) ? remake(r; nominal_rate=r.nominal_rate*k) : r for (i, r) in enumerate(m.reactions)])
scale_kpred(k)  = m -> remake_model(m; reactions=Reaction[
    i==7 ? remake(r; half_saturation_constants=Dict("HET"=>r.half_saturation_constants["HET"]*k)) : r
    for (i, r) in enumerate(m.reactions)])
set_water_factor(v) = m -> remake_model(m; reactions=Reaction[
    i==7 ? remake(r; efficiency_flowing=v) : r for (i, r) in enumerate(m.reactions)])

# filter perturbation builders (f -> f', filter is mutable)
scale_light_sand(k)  = f -> (f.light_attenuation_sand *= k; f)
scale_light_water(k) = f -> (f.light_attenuation_water *= k; f)
set_temp(v)          = f -> (f.temperature = v; f)

# ---- headline output metrics ------------------------------------------------
function metrics(r)
    z = depths(r); dz = gridsize(r.filter); poro = computeporosity(r.filter, z)
    si = findfirst(>=(0.0), z)                    # first sand cell
    matHET = concentration(r, "HET", :matrix)[:, end]
    flowPAT = concentration(r, "PAT", :flowing)[:, end]
    o2 = concentration(r, "O2", :flowing)[:, end]
    attHET = sum(poro .* matHET .* dz)
    attPAT = sum(poro .* concentration(r, "PAT", :matrix)[:, end] .* dz)
    surf = matHET[si]
    pen_i = surf > 0 ? findlast(>(0.01*surf), matHET) : nothing
    pen = pen_i === nothing ? 0.0 : z[pen_i]
    eff = flowPAT[end]
    removal = 1 - eff / INFLOW[4]
    o2_sup = o2[max(1, si-3)]                      # O2 just above the sand
    maxphib = maximum(get_volume_fractions(r).biofilm[:, end])
    return (; flag=r.flag, attHET, attPAT, pen, eff, removal, o2_sup, maxphib)
end

function run_case(fmod, mmod)
    f = fmod(base_filter())
    m = mmod(base_model())
    s = State(f, m); seed_mature!(s, f.grid.centers)
    r = simulate(s; inflow_concentrations=INFLOW, simulation_time=TSIM,
                 time_step=DT, n_frames=NFRAMES, clogging_fraction=0.99, quiet=true)
    return metrics(r)
end

# ---- perturbation list: (label, group, filter-mod, model-mod) ----------------
id = x -> x
const PERTURB = [
    # attachment / detachment
    ("attachment ×0.5",   "attach/detach",  id, scale_attach(0.5)),
    ("attachment ×2",     "attach/detach",  id, scale_attach(2.0)),
    ("sand_pathogen=0.1", "attach/detach",  id, set_sand_pat(0.1)),
    ("sand_pathogen=0.5", "attach/detach",  id, set_sand_pat(0.5)),
    ("detachment ×0.5",   "attach/detach",  id, scale_detach(0.5)),
    ("detachment ×2",     "attach/detach",  id, scale_detach(2.0)),
    # biofilm excess velocity (Cahn-Hilliard)
    ("zeta_0 ×0.5",       "biofilm vel.",   id, scale_zeta0(0.5)),
    ("zeta_0 ×2",         "biofilm vel.",   id, scale_zeta0(2.0)),
    # light attenuation
    ("light_sand ×0.5",   "light",          scale_light_sand(0.5),  id),
    ("light_sand ×2",     "light",          scale_light_sand(2.0),  id),
    ("light_water ×0.5",  "light",          scale_light_water(0.5), id),
    ("light_water ×2",    "light",          scale_light_water(2.0), id),
    # environmental reaction rates
    ("growth ×0.5",       "reaction rate",  id, scale_rx((1,2), 0.5)),
    ("growth ×2",         "reaction rate",  id, scale_rx((1,2), 2.0)),
    ("all rates ×0.5",    "reaction rate",  id, scale_rx((1,2,3,4,5,6,7), 0.5)),
    ("all rates ×2",      "reaction rate",  id, scale_rx((1,2,3,4,5,6,7), 2.0)),
    ("temperature=10°C",  "reaction rate",  set_temp(10.0), id),
    ("temperature=30°C",  "reaction rate",  set_temp(30.0), id),
    # PAT inactivation / bacterivory
    ("inactivation ×0.5", "PAT kinetics",   id, scale_rx((6,), 0.5)),
    ("inactivation ×2",   "PAT kinetics",   id, scale_rx((6,), 2.0)),
    ("bacterivory ×0.5",  "PAT kinetics",   id, scale_rx((7,), 0.5)),
    ("bacterivory ×2",    "PAT kinetics",   id, scale_rx((7,), 2.0)),
    ("kPred ×0.5",        "PAT kinetics",   id, scale_kpred(0.5)),
    ("kPred ×2",          "PAT kinetics",   id, scale_kpred(2.0)),
    ("water_factor ×10",  "PAT kinetics",   id, set_water_factor(1e-2)),
    ("water_factor ×0.1", "PAT kinetics",   id, set_water_factor(1e-4)),
]

# ---- drive ------------------------------------------------------------------
println("Base case: mature-seeded pathogen filter, T=$TSIM d, dt=$DT, $NCELL intervals")
b = run_case(id, id)
@printf("Base metrics: flag=%s  attachedHET=%.4e  attachedPAT=%.4e  PATeffluent=%.4e  O2above=%.4e  maxφb=%.3f\n\n",
        b.flag, b.attHET, b.attPAT, b.eff, b.o2_sup, b.maxphib)

reldiff(a, base) = base == 0 ? (a == 0 ? 0.0 : Inf) : (a - base) / abs(base)

rows = Vector{Any}[]
push!(rows, ["label","group","flag","attachedHET","attPAT","PATeffluent","O2above","maxphib",
             "d_attHET","d_attPAT","d_PATeff","d_O2"])
sens = Tuple{String,String,Float64,Float64,Float64,Float64}[]  # label, group, |Δ| for each key metric
for (label, group, fmod, mmod) in PERTURB
    mt = run_case(fmod, mmod)
    d_att  = reldiff(mt.attHET, b.attHET)
    d_pat  = reldiff(mt.attPAT, b.attPAT)
    d_eff  = reldiff(mt.eff,    b.eff)
    d_o2   = reldiff(mt.o2_sup, b.o2_sup)
    push!(rows, [label, group, mt.flag, mt.attHET, mt.attPAT, mt.eff, mt.o2_sup, mt.maxphib,
                 d_att, d_pat, d_eff, d_o2])
    push!(sens, (label, group, d_att, d_pat, d_eff, d_o2))
    @printf("  %-20s [%-13s] flag=%-8s  ΔattHET=%+.1f%%  ΔattPAT=%+.1f%%  ΔPATeff=%+.1f%%  ΔO2=%+.1f%%\n",
            label, group, mt.flag, 100d_att, 100d_pat, 100d_eff, 100d_o2)
end

# write CSV
outcsv = joinpath(@__DIR__, "sensitivity_results.csv")
writedlm(outcsv, rows, ',')
println("\nwrote $outcsv")

# ranked summary: mean |relative change| across the four key metrics
println("\nMost influential perturbations (mean |Δ| over attachedHET, attPAT, PATeffluent, O2):")
score(s) = (m = filter(isfinite, [abs(s[3]), abs(s[4]), abs(s[5]), abs(s[6])]); isempty(m) ? 0.0 : sum(m)/length(m))
for (label, group, d1, d2, d3, d4) in sort(sens; by=score, rev=true)[1:min(10, length(sens))]
    @printf("  %6.1f%%  %-20s (%s)\n", 100*score((label,group,d1,d2,d3,d4)), label, group)
end
