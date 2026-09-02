# Candidate figures for AWR-SSF

Plots that exist in the code and could earn a place in the manuscript, but are not in it.
Each entry: what it shows, how to make it, what it would argue, and what has to be settled
before it could be used. Nothing here is committed to the paper.

---

## 1. Biofilm composition against depth (`Results.plotBiofilmComposition`)

**What it shows.** A stacked-area plot of the biofilm's *relative* composition at one
instant: what fraction of the biofilm at each depth is HET, PHO, POM, PAT, by volume
fraction (mass / component density). `Phase="Biofilm"` sums matrix + enclosed;
`Phase="Matrix"` uses the matrix alone.

**How to make it.**

```matlab
D = load('analysis/probes/data/chain/chain_fld2x_lit_leg3.mat','results');
plotBiofilmComposition(D.results, 30, FontSize=13);   % t = 30 d
```

**What it would argue.** The paper currently reports total biomass and per-component
profiles separately; nothing shows the *community structure* varying with depth. On the
E4 lit chain at 30 d the plot reads:

- HET is 62-88% of biofilm volume throughout, rising monotonically with depth;
- PHO is 12-37%, peaking near the surface (maximum around z = 0.04 m, not at z = 0);
- POM and PAT are slivers at this scale, essentially invisible.

That is a one-panel statement of the Schmutzdecke argument the discussion makes in prose,
and it complements `fig:pat-scraping`: it shows *why* removing the top layer matters less
than expected, because heterotrophs — which drive bacterivory — dominate at every depth.

**Before it could be used.**

- The PHO maximum sits *below* the surface. Worth understanding before publishing: light
  is greatest at z = 0, so a subsurface peak needs an explanation (self-shading?
  attachment? the roughness layer?), not just a caption.
- PHO stays ~11% at z = 1 m, where no light reaches and `MinimumLightFactor = 0` since
  2026-08-26 — so phototrophs there cannot be growing. That fraction must arrive by
  transport/attachment. Fine, but the caption must not imply growth.
- Relative fractions hide absolute magnitude: at 30 d the whole column holds
  0.221 kg/m^2 of particulate biomass. Pair the plot with an absolute panel, or state
  phi_b in the caption. (CORRECTED 2026-09-02: this entry first said "phi_b max 0.0032 ...
  a nearly-empty filter". That was the PARTICULATE volume fraction mislabelled as phi_b.
  True phi_b = 0.311 at z = 0, 0.145 at 1 cm, 2.2e-3 at 0.5 m, 2.3e-5 at 1 m -- exactly
  100x larger, because biofilm is 99 % enclosed water at beta = 0.99. The schmutzdecke is
  NOT nearly empty; the deep bed is. The pairing advice stands, the number did not.)
- Would need re-running under the Tenore presets (b1fcd03) like everything else.

**Status.** Method fixed 2026-09-02 (see below). Not yet a manuscript figure.

---

## Method fixes applied 2026-09-02

`src/@Results/plotBiofilmComposition.m`, written 2026-09-02 17:52, had four defects:

1. **No `hold`** — `area()` clears the axes by default, so the loop drew four bands and
   kept only the last. The plot showed ONE component, not a stack. This is the one that
   mattered; the others are cosmetic.
2. **`legend()` with no handle** targets `gca`, which is not `ax` when the caller passes
   `AxisHandle` — the legend landed on whatever axes happened to be current.
3. **`FontSize` and `Grid` declared but never applied.**
4. **`Phase` had no `mustBeMember`**, so a typo fell through the `switch` and left
   `particle_biofilm_mass` undefined, erroring later at a confusing place.

Verified after the fix: 4 area objects drawn, legend entries PAT/POM/PHO/HET parented to
the supplied axes, FontSize and grid applied, and `Phase="Flowing"` now rejected at the
argument block. Axis labels added.

---

## 2. Community structure against depth, four scenarios (`recreateFigsManriquez2026("composition")`)

**What it shows.** Eight panels, 2 x 4. Top row: the stacked-area composition of
`plotBiofilmComposition` (relative volume fraction HET/PHO/POM/PAT, matrix + enclosed)
against depth over the whole bed 0 <= z <= 1 m, one panel per scenario — Summer
(`fld2x_lit`, 19 C uncovered), Winter (`fld2x_winter`, 3 C), Covered (`fld2x_cov01`,
1 % light), Dark (`fld2x_dark`, 0 % light) — all at the *same* instant. Bottom row: the
absolute measure the top row hides, phi_b(z) on a log axis, same depth axis, common
y-limits across scenarios.

The horizon is chosen, not assumed: `scenarioSet` takes the latest time ALL four chains
reach. With the chains on disk on 2026-09-02 that is **t = 30 d** (`fld2x_cov01` stops at
leg 3; lit and dark run to 104 d, winter to 90 d).

**How to make it.**

```matlab
cd analysis
recreateFigsManriquez2026("composition")   % -> results/figures/recreation/rec_composition_depth.{png,pdf}
```

**What it would argue.** Two things at once, which is why it is worth a whole figure.

1. *Community structure is a depth structure, and it barely responds to light.* HET rises
   monotonically from ~62 % at the surface to ~89 % at z = 1 m in summer; PHO is the
   complement. Summer, Covered and Dark are visually indistinguishable. Winter is the
   only panel that differs: HET drops to 34 % at z = 0 and PHO takes 66 % of the
   schmutzdecke.
2. *The relative picture must not be read alone.* The bottom row shows phi_b falling four
   decades across the bed (0.311 at z = 0, 0.145 at 1 cm, 2.2e-3 at 0.5 m, 2.3e-5 at
   1 m, summer at 30 d). The 89 % HET at the bottom of the bed describes essentially
   nothing.

Together that is the Schmutzdecke argument in one figure, and it complements
`fig:pat-scraping`: what scraping removes is the only part of the column with biomass in
it, but the community there is not compositionally special except in winter.

**What must be settled first.**

- **The chains are stale.** All four `fld2x_*` chains were grown under the OLD presets
  (mu_HET = 2.0, nu_P = 0.094). The current presets (b1fcd03) carry mu_HET = 4.8 and
  nu_P = 52. Nothing here can go in the paper until the four chains are regenerated
  under b1fcd03 — and mu_HET more than doubling is not a perturbation.
- **`fld2x_cov01` is short.** It stops at 30 d, which is what caps the whole figure. The
  paper discusses 90 d behaviour; either extend cov01 or state the 30 d horizon plainly.
- **PHO at depth is not growth.** I/I0 < 1e-5 below 1 cm (sand attenuation), and
  `MinimumLightFactor = 0` since 2026-08-26, so the phototrophs at z = 0.5 m cannot be
  photosynthesising. They arrive by transport and attachment. The caption must say so.
- The subsurface PHO maximum noted in entry 1 is still unexplained.

**Status.** Built 2026-09-02 (agent B, autonomous block). Not a manuscript figure.

---

## 3. Bed-integrated composition over time, Tenore2021 style (`recreateFigsManriquez2026("tenore")`)

**What it shows.** Four panels, one per scenario. In each, a stacked area of the
*fractions* f_i(t) of sessile (matrix + enclosed) biomass integrated over the bed
0 <= z <= 1 m — HET, PHO, POM, PAT — from t = 0 to the common horizon, with the
**absolute** areal biomass sum_i m_i(t) drawn as a black line on the right axis, in
kg/m^2 WET (concentrations as the model carries them, rho_P = 1117).

**Why this presentation, and why it is the Tenore analogue.** Tenore, Mattei and Frunzo
(2021, *Commun Nonlinear Sci Numer Simulat* 94:105577) present phototroph-heterotroph
competition as *fractions f_i evolving in time* — their Fig. 2 is biofilm thickness L(t),
Fig. 3 the per-species sessile mass at T = 10, 15, 20 d, Fig. 6 the phototrophic mass at
four times. Their whole argument (heterotrophic pioneers facilitate phototrophic
invasion) is a statement about how the *mix* moves, not about a profile at one instant.
Our model has no free biofilm surface, so their L(t) has no counterpart here — the
1-D continuum filter has no thickness to plot. What it does have is the same competition,
resolved over the bed, and the honest analogue is the bed-integrated f_i(t). It is the
time-analogue of entry 2 and the manuscript has nothing like it: every existing figure is
either a profile at one instant or a scalar time series (effluent, total mass). None
shows the community *re-sorting* itself.

**How to make it.**

```matlab
cd analysis
recreateFigsManriquez2026("tenore")   % -> results/figures/recreation/rec_composition_time.{png,pdf}
```

**What it would argue.** At t = 30 d, bed-integrated:

| scenario | areal biomass [kg/m^2] | HET | PHO | POM | PAT |
|---|---|---|---|---|---|
| Summer  (`fld2x_lit`)    | 0.0845 | 0.653 | 0.335 | 0.012 | 0.000 |
| Winter  (`fld2x_winter`) | 0.1924 | 0.543 | 0.447 | 0.010 | 0.000 |
| Covered (`fld2x_cov01`)  | 0.0825 | 0.656 | 0.333 | 0.012 | 0.000 |
| Dark    (`fld2x_dark`)   | 0.0825 | 0.656 | 0.333 | 0.012 | 0.000 |

- **Covered and dark are the same run to four significant figures.** Between 0 % and 1 %
  incident light nothing in the bed changes. That is a real, publishable negative: the
  light pathway does not reach the bed, and any covered/uncovered claim the paper makes
  has to be about the supernatant and the schmutzdecke surface, not the filter bed.
- **PHO is a third of bed biofilm in total darkness.** Same point as entry 2, but here it
  is unarguable — the dark panel has no light at all and PHO still stabilises at 0.333.
  Phototrophs in the bed are an *attachment* signal from the influent.
- **Winter accumulates 2.3x the summer biomass and is still climbing** at 30 d while
  summer has plateaued (detachment balancing growth). Consistent with Bae2023 (winter
  >= summer) noted 2026-08-25.
- The composition is essentially settled by ~t = 15 d in every arm; the transient is the
  first two weeks.

**What must be settled first.**

- **Stale presets, exactly as entry 2** — mu_HET 2.0 -> 4.8 and nu_P 0.094 -> 52 (b1fcd03).
  The HET/PHO split is precisely the quantity a 2.4x change in mu_HET will move, so this
  figure is the *most* preset-sensitive of the three candidates. Regenerate before use.
- **The 30 d cap comes from `fld2x_cov01`.** If cov01 is not extended, either drop the
  covered arm and run to 90 d, or keep four arms at 30 d. Do not silently mix horizons.
- **Wet mass.** The right axis is wet areal biomass. If the paper reports dry mass
  elsewhere, one of the two has to be converted; state the convention in the caption.
- **PAT is identically ~0** at this influent, so the fourth band is invisible. Either drop
  PAT from the stack or say in the caption why it is absent.
- The covered/dark identity is physics, **not** a preset collision: `rec.lightScale` reads
  1 (lit), 0.01 (cov01), 0 (dark) in leg 3 of each chain, and phi_b at z = 0 differs
  (0.3114 / 0.2708 / 0.2701). Checked 2026-09-02. Note the wider point that follows:
  full light versus 1 % light changes the bed-integrated areal biomass only from 0.0845
  to 0.0825 kg/m^2 — a 2.4 % effect for a 100x change in light. The bed is dark whatever
  the cover does; the light-dependent action is all in the supernatant and the top
  millimetres.

**Status.** Built 2026-09-02 (agent B, autonomous block). Not a manuscript figure.
