# Recreation of the Manriquez2026 results figures

Written by `analysis/recreateFigsManriquez2026.m`. Regenerate everything with

```
cd /Users/jaime/Research/SSF/code/1d/MPC-SSF
/Applications/MATLAB_R2025b.app/bin/matlab -batch "addpath('analysis'); recreateFigsManriquez2026('all')"
```

or one family at a time: `recreateFigsManriquez2026("seasons")`, keys
`light | seasons | roofed | month | scrape | outflow1d | outflow2d | patpulse |
filtration | pat2d`.

**What this is.** Each panel copies the published figure's *design* — panel
structure, axis orientation, ranges, tick spacing, labels, units, legend text
and placement, line styles, marker sets, zoom windows — and fills it with the
audited working-set data of this repo (E4/E7: 2x field influent, zeta0 = 1,
zeta1 = 0.27, linear detachment, transfer x10, N = 500, 19 C). **The curves are
not expected to agree with the published ones**; only the framing is meant to
be one-to-one comparable. Our biofilm is roughly an order of magnitude lighter
than the published one (0.089 vs 1.19 kg/m^2 at steady state) and peaks at
phi_b ~ 0.33 rather than ~ 0.56, so several published axis windows had to be
slid; every such slide is listed below and keeps the published window *width*.

Published figures live in
`/Users/jaime/Research/SSF/manuscripts/AWR-SSF/figures/`.

## Figure table

| key | published file | recreated file | data source | deliberate deviation |
|---|---|---|---|---|
| light | `results_light-seasons-time.pdf` | `rec_light-seasons-time` | analytic | third dotted curve "Working set" added: the light curve actually wired into probeChain, so the reader can see how our forcing sits between the two published ones. Published Summer/Winter curves reproduced exactly. |
| seasons | `results_biofilm_seasons.pdf` | `rec_biofilm_seasons` | `chain_fld2x_lit_leg1..9` at t = 90 d | Complete since the 2026-09-01 re-render: both arms at t = 90 d (`chain_fld2x_winter_leg1..9` arrived from cosmos job 3550889 on 08-28). |
| seasons | `results_biofilm_seasons_zoomed.pdf` | `rec_biofilm_seasons_zoomed` | same, z in [-0.02, 0.04] | same Winter omission. |
| roofed | `results_biofilm_covered.pdf` | `rec_biofilm_covered` | `chain_fld2x_dark` (Covered) vs `chain_fld2x_lit` (Uncovered), both at t = 30 d | Complete since the 2026-09-01 re-render: "Covered" is the 1 % arm (`fld2x_cov01`, arrived from cosmos job 3550892 on 08-28), matching the published 0.1 %-of-light setup more closely than the 0 % fallback. |
| roofed | `results_biofilm_covered_zoomed.pdf` | `rec_biofilm_covered_zoomed` | same, z in [-0.025, 0.01] | same covered-arm fallback. |
| roofed | `results_biofilm_covered-ultrazoomed_zoomed.pdf` | `rec_biofilm_covered-ultrazoomed_zoomed` | same, z in [0.10, 0.11] | phi window re-centred on our own curves. The published window is phi in [0.39, 0.405]; our phi_b at that depth is far lower, so a fixed window would show empty axes. The window **width** (0.015) and the depth window ([0.10, 0.11]) are the published ones; only the centre moves. Same covered-arm fallback. |
| month | `results_evolution_summer.pdf` | `rec_evolution_summer` | `chain_fld2x_lit_leg1..2`, t = 1, 2, 5, 10, 20 d | none. |
| month | `results_evolution_summer_zoomed.pdf` | `rec_evolution_summer_zoomed` | same, z in [-0.02, 0.04] | none. |
| month | `results_evolution_summer_mass.pdf` | `rec_evolution_summer_mass` | `chain_fld2x_lit_leg1..9`, 0-90 d | published y range is [0, 1.2] kg/m^2; ours is autoscaled to [0, 1.05 x max] because our steady-state mass is 0.089 kg/m^2 and the published range would flatten the curve onto the axis. x range 0-90 d is the published one. |
| scrape | `results_evolution_Scraped_GP0_zoomed.pdf` | `rec_evolution_Scraped_GP0_zoomed` | `analysis/results/figures/repro/data/scrape_GP0.mat` | Complete since the 2026-09-01 re-render on the 20 d reruns (cosmos job 3550884, on disk 08-28); the t = 50 d curve is restored. |
| scrape | `results_evolution_Scraped_GP4_zoomed.pdf` | `rec_evolution_Scraped_GP4_zoomed` | `scrape_GP4.mat` | as above. |
| scrape | `results_evolution_Scraped_GP8_zoomed.pdf` | `rec_evolution_Scraped_GP8_zoomed` | `scrape_GP8.mat` | as above. |
| scrape | `results_evolution_Scraped_GP12_zoomed.pdf` | `rec_evolution_Scraped_GP12_zoomed` | `scrape_GP12.mat` | as above. |
| scrape | `results_evolution_Scraped_GP_mass.pdf` | `rec_evolution_Scraped_GP_mass` | `scrape_GP{0,4,8,12}.mat` | x axis kept at the published 30-60 d, so the right half is empty until the 20 d reruns land; y range autoscaled (published 0.4-1.2 kg/m^2, ours 0.03-0.09). |
| outflow1d | `Outflow1D{O2,IC,NH4,HPO4,DOM}.pdf` | `rec_Outflow1D{O2,IC,NH4,HPO4,DOM}` | `chain_fld2x_lit_leg1..3`, effluent rows, t <= 30 d | y range autoscaled per species (published axes are per-species too); records are kg/m^3 = g/L and are multiplied by 1000 for mg/L. |
| outflow2d | `Outflow2D{O2,IC,NH4,HPO4,DOM}.pdf` | `rec_Outflow2D{O2,IC,NH4,HPO4,DOM}` | `chain_fld2x_lit_leg1..3`, `results_py.liquids(:,:,j,2)` (flowing suspension), t <= 30 d | z (concentration) range autoscaled; axes, view, jet colormap and colorbar label match the published panels. |
| patpulse | `March-11/PATPulse_HighLow.pdf` | `rec_PATPulse_HighLow` | `pulse_fig_p1x_d37.mat` (Pulse 1), `pulse_fig_p1em3_d37.mat` (Pulse 2) | none. IN curves are reconstructed as `cPeak` over `[PulseT0, PulseT1]` (the record stores only the effluent). |
| patpulse | `March-11/PATPulse_FastSlow.pdf` | `rec_PATPulse_FastSlow` | `pulse_fig_{p1x,exp2,exp3}_d37.mat` | none. |
| filtration | `March-11/FiltrationRate_HETPHO.pdf` | `rec_FiltrationRate_HETPHO` | `pulse_fig_{p1x,p1em3}_d37.mat` | y window slid from the published [4, 7] to cover our removals (HET ~ 3.0, PHO ~ 4.5); the published span of three decades and the 0.5 tick spacing are kept. Like the published panel, the axes carry no labels. |
| pat2d | `March-11/PAT{Flowing,Matrix,Enclosed}Reference.pdf` | `rec_PAT{Flowing,Matrix,Enclosed}Reference` | `pulse_fig_p1x_d37.mat`, `results.Frames.Concentrations{"PAT", vol}` | colour limits autoscaled per panel (log scale, jet, as published). Skipped with a printed line if the file holds no `results` object. |

## One trap worth naming

`rec.effluent` rows follow `rec.effNames`
(`O2 IC NH4 HPO4 DOM HET PHO POM PAT`), but `rec.influent` / `rec.Influent`
follows probeChain's `InflowConcentrations` argument order
(`HET PHO POM PAT O2 IC NH4 HPO4 DOM`). Indexing the influent with the
effNames position gives HET = 6.23e-3 instead of 3.0e-4 and inflates the HET
log removal by ~1.3 decades. `influentIndex` in the script exists for exactly
this; `liquidIndex` in `reproduceManriquez2026.m` uses the same convention.

## Relation to `../repro/`

`reproduceManriquez2026.m` draws the same experiments in a *diagnostic* style
(tiled multi-panel figures, colour, extra series). This script is the
complementary one: one file per published panel, in the published style, for
side-by-side visual comparison. Neither replaces the other and neither touches
the other's outputs.

## `rec_BiologicalActivity` — activity from the reaction terms (not a published panel)

`analysis/plotBiologicalActivity.m` writes `rec_BiologicalActivity.{png,pdf}`
(add `Split=true` for `rec_BiologicalActivity_profiles` and
`_timeseries`). It has no published counterpart: it measures *biological
activity* directly from the model's reaction terms rather than from a standing
stock. For each reaction *j* it evaluates the volumetric extent rate
`r_j = phi * mu_j(T) * I_j(z) * min_k Monod_k * X_driver_local` in each of the
three phases the solver uses, by copying the assembly of `src/@State/simulate.m`
(lines 70-84 `listK`/`listOrder`, 103-137 porosity/`dz`/eta, 374-401 `muRates`,
`lightOptimal`, light-mode flags and phase efficiencies, 415-462 global -> local
concentrations, 472-497 attenuation and `lightFactor`, 506-514 `ecoRx*`, and the
local `evaluateReactions` at 1115) onto the `snap` state of each
`chain_fld2x_{lit,dark}_leg{1..11}.mat`. The model and filter are rebuilt with the
E4/E9 working set exactly as `probeChain.m:74-120` does (`Zeta0=1`, `Zeta1=0.27`,
`KDOM=3e-4`, `KHPO4=1e-6`, `TransferScale=10`, linear detachment, `EtaSand=1500`,
`Delta=5e-3`, 19 C, `LightForm="chain"`). Panel (a) is the depth profile of the
rate at 104 d over z in [-0.02, 0.30] m, panel (b) the column-integrated rate
`int eps(z) r_j dz` at each leg time 10..104 d — the `eps` weight is the solver's,
since `simulate.m:966,984` advances a quantity whose conserved form is `eps*g`.

**Assumptions, stated because they are load-bearing.** (i) `snap.Time` is a whole
number of days, where the chain light curve is exactly zero, so a snapshot-instant
evaluation would report zero phototroph growth. The state is therefore frozen and
the rates are evaluated at 48 times spanning one day; every number plotted is a
**daily mean**. This is only legitimate because the state is quasi-steady by 30 d
(E9: 0.001 %/10 d profile change) — it would be wrong during ripening, and panel
(b)'s 10-30 d points inherit that caveat. (ii) The PAT-driven reactions
(MarkerGrowth, Inactivation, Bacterivory) are identically zero because PAT is
absent from the E4/E9 influent, and are omitted. (iii) The flowing phase
contributes nothing: `EfficiencyFlowing = 0` for every reaction except
Bacterivory (`simulate.m:394-396`), so 99.997 % of the column activity is
biofilm-matrix and 3e-5 is the enclosed phase. Suspension in the supernatant is
inert by construction, which is why every profile is zero above z ~ -6 mm.

**Validation.** Column-integrated net O2 reaction sink at 104 d, lit:
**21.115 g O2/(m^2 d)** against the advective deficit `q(C_in - C_out)` =
**21.852 g O2/(m^2 d)** (`q` = 7.2 m/d = `SandFilter.InflowVelocity`, C_in =
9.10e-3, C_out = 6.065e-3 kg/m^3), a **-3.4 %** difference — within the margin
left by the diel cycle (the sink is a daily mean, C_out a midnight value) and by
outlet dispersion, which the advective estimate ignores. Phototroph growth is
**exactly 0** everywhere in the dark arm, and in the lit arm is confined to
z in [-0.004, 0.004] m. Steady-state column activity, g/(m^2 d), lit vs dark:
hydrolysis 24.40 / 23.60, HET death 21.56 / 21.05, HET growth 19.41 / 18.89,
PHO death 7.73 / 7.20, PHO growth 0.794 / 0. So light adds 3-7 % to every
heterotrophic term and is the sole source of the 0.79 g/(m^2 d) of phototroph
production; the filter's activity is overwhelmingly heterotrophic even in the
lit arm.

## `rec_supernatant_height.{pdf,png}` — biofilm height in the supernatant

Not a recreation of a published panel: there is no published figure of this. Added
2026-09-02 to put a number on the `results.tex` caption claim that winter and summer
have "comparable biofilm heights above" the bed.

**Measure.** `analysis/plotSupernatantHeightChains.m`, backed by
`src/@Results/getSupernatantHeight.m`:

    h(t) = ( \int_{z<0} eps(z) phi_b(z,t) dz ) / phi_ref,      phi_ref = 0.3

the thickness the biofilm standing above the sand surface (z = 0) would have if
compacted to volume fraction `phi_ref`. `phi_ref = 0.3` is the biofilm volume
fraction at the sand surface in the working set (0.271 dark, 0.271 cov01, 0.306
winter, 0.324 summer at their final frames). `h` is linear in `1/phi_ref`, so the
ordering of the scenarios does not depend on that choice.

**Why not a threshold height.** The obvious measure — the most negative z with
phi_b above a threshold — cannot be used at N = 500. The supernatant biofilm
occupies one or two cells (dz = 1.998 mm), so the threshold height can only return
0, 2 or 4 mm; at thr = 1e-3 all four scenarios read 4 mm, and winter drops from
4 mm to 2 mm on the same frame when the threshold moves to 1e-2.

**Caveat.** h is of order 0.5-0.9 mm, i.e. h/dz = 0.27-0.43. The supernatant
biofilm is *not resolved* at N = 500 and every number here is a sub-cell quantity.
The existing chain diagnostics `rec.supInt` / `rec.supFrac` are identically **zero**
on these chains, because they take the supernatant to be z < -SandRoughness =
z < -5 mm and all of the biofilm above the sand sits inside the roughness band.

**Numbers** (final frame; summer/dark run to 104 d, winter to 90 d, cov01 to 30 d):

| scenario | h(30 d) | h(90 d) | h(final) | bed integral at 90 d |
|---|---|---|---|---|
| summer, `fld2x_lit`    | 0.778 mm | 0.861 mm | 0.861 mm | 0.00772 |
| winter, `fld2x_winter` | 0.636 mm | 0.679 mm | 0.679 mm | 0.02411 |
| covered 1 %, `fld2x_cov01` | 0.538 mm | — | 0.538 mm | — |
| covered dark, `fld2x_dark`  | 0.534 mm | 0.541 mm | 0.541 mm | — |

`cov01` and `dark` overlie each other to within 1 % and are hard to tell apart on
the figure; 1 % light is indistinguishable from none by this measure.

**Reading.** Summer is the *taller* arm: 0.861 vs 0.679 mm at 90 d, +27 %. Winter
carries 3.12x the bed biomass (0.02411 vs 0.00772). Winter leads summer in height
only over the first ~15 d. So the caption's "comparable" is fair on magnitude
(27 % apart, against 212 % apart in the bed) but the *sign reverses* between bed
and supernatant, and the taller arm is summer.
