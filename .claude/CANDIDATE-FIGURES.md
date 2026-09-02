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
- Relative fractions hide absolute magnitude: at 30 d the whole column holds only
  0.221 kg/m^2 (phi_b max 0.0032). A composition plot of a nearly-empty filter can
  mislead. Pair it with a total-biomass panel, or state phi_b in the caption.
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
