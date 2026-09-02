# Sand-attachment sweep with scraping: can we reproduce Schijven2013's scraping loss?

## Goal

Schijven2013 (experiments D1 -> D2) scraped the Schmutzdecke off a ripened slow sand filter
and measured **0.6 log10 lower removal of MS2 and 1.6 log10 lower removal of E. coli WR1**.
A 53 d Schmutzdecke gave 0.8 and 2.2 log10 higher removal than a young one.

Our marker currently has SandAttachmentFactor s = 0: it cannot attach to bare sand, only to
biofilm. That is the extreme of "biofilm is the only collector". Schijven's own model has a
sand-attachment baseline PLUS a Schmutzdecke term, so s = 0 removes the baseline he fits.

Done means: log removal measured on an intact and a freshly scraped 20 d filter, for a
sweep of s, and the scraping loss `Delta L = L_intact - L_scraped` compared against
Schijven's 0.6-1.6 log.

**The prediction to test.** At s = 0 all removal is biofilm-mediated, so scraping should cost
the MOST. As s rises the bare sand carries more of the removal and scraping should cost LESS.
If the measured Delta L brackets 0.6-1.6 log somewhere in the sweep, that value of s is the
one consistent with Schijven, and it is evidence on the question of whether sand or biofilm
attachment should be larger.

## Sweep

s in {0, 0.06, 0.3, 1.0}. 0 is the current model; 0.06 is Schijven2013's T4 sticking
efficiency (geometric mean), already the OAT nominal; 1.0 is "sand as sticky as biofilm",
the Diehl2025 formulation with no asymmetry; 0.3 spans the middle.

Base filter: `chain_fld2x_lit_leg2`, the 20 d working-set filter (T = 20 d as asked).
Scrape depth 2 cm, the depth `fig:pat-scraping` uses and the one carrying 18 % of removal
activity (measured 2026-09-02).

## Steps

1. `analysis/probes/probePulse.m` — add `SandPathogen` option (default NaN = preset 0),
   forwarded to `pathogenModel`. Verify: `probePulse(..., SandPathogen=0.06)` gives PAT
   SandAttachmentFactor 0.06 in the saved `results.Model`.
2. `analysis/reproduceManriquez2026.m` — allow `Days = 0` in `runScrape`: scrape the state
   and save the snapshot WITHOUT simulating, so a filter can be challenged immediately
   after scraping as in Schijven's D2. Verify: a `Days=0` call writes a chain snapshot whose
   biofilm above the scrape depth is zero and whose `Time` equals the parent leg's.
3. `slurm/sandfactor_scrape.sbatch` — array 0-7: {intact, scraped} x s. Each arm is one
   `probePulse` on the appropriate snapshot at `TPost = 3` d.
4. Compare `rec.Lmin` / the removal curves across arms; tabulate `Delta L` per s.

## Files touched

- `analysis/probes/probePulse.m`, `analysis/reproduceManriquez2026.m`,
  `slurm/sandfactor_scrape.sbatch`. No `src/` edit: `pathogenModel` already exposes
  `SandPathogen`.

## Caveat to carry into any conclusion

Schijven measured MS2 and E. coli; our marker is neither, and its inactivation and
bacterivory rates are Manriquez Table B.4 values, not fitted to his organisms. A match in
Delta L would be evidence about the ATTACHMENT structure, not a validation of the whole
pathogen submodel.
