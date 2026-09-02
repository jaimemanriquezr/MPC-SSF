# Plan — a measure of biofilm height in the supernatant

Autonomous block 2026-09-02, agent A. Brief: `.claude/autonomous/2026-09-02-newparams.md`
(workspace `.claude/`).

## Goal

Quantify how far biofilm reaches into the supernatant water (z < 0, above the sand
surface) and plot it, so the manuscript claim "comparable biofilm heights above it"
can be tested against the working set instead of asserted.

Done means: a getter and a plot method on `@Results`, a figure for the four `fld2x_*`
chains in `analysis/results/figures/recreation/`, and the final heights reported.

## The measure

Two candidates were evaluated on the committed chains (`rec.phiT`, which is exactly
`getVolumeFractions().Biofilm`):

(a) **Threshold height** — the most negative z with phi_b > thr.
    REJECTED. At N = 500 (dz = 1.998 mm) the supernatant biofilm occupies one or two
    cells, so h_thr can only take the values 0, 2 or 4 mm. It cannot separate the
    scenarios (summer, winter, dark and cov01 all read 4 mm at thr = 1e-3) and it
    flips on the threshold (winter: 4 mm at thr = 1e-3, 2 mm at thr = 1e-2, same frame).

(b) **Equivalent height** — CHOSEN.

        h_eq(t) = ( \int_{z<0} phi_b(z,t) dz ) / phi_ref

    the thickness the supernatant biofilm would have if compacted to a reference
    volume fraction. Default `phi_ref = 0.3`, the biofilm volume fraction at the sand
    surface in the working set (0.271 dark / 0.271 cov01 / 0.306 winter / 0.324 summer
    at the final frame). No threshold; continuous in time; the ranking between
    scenarios is independent of phi_ref because h_eq is a pure rescale of the
    supernatant integral (`rec.supInt` with porosity, which is 1 above the bed).

Porosity is not applied: above z = 0 there is no sand, eps == 1 there (verified on
`rec.eps`), so `sum(phi_b) * dz` and `sum(eps .* phi_b) * dz` agree in the supernatant.
The method still multiplies by porosity so that a non-unit-porosity supernatant, or a
caller passing a bed window, stays correct.

## Steps

1. `src/@Results/getSupernatantHeight.m` — `[h, t] = getSupernatantHeight(obj, options)`
   with `options.Reference (1,1) double = 0.3` and `options.SurfaceDepth (1,1) double = 0`
   (the z that separates supernatant from bed). Returns h [1 x nT] in metres and the
   frame times.
   *Verify:* on `chain_fld2x_lit_leg9.mat` the final value matches the Python reading
   of `rec.supInt / 0.3` to 4 s.f. (1.2865 mm).
2. `src/@Results/plotSupernatantHeight.m` — time series, styled like
   `plotBiofilmComposition.m` (AxisHandle / FontSize / Grid options that are all
   actually used, `hold` guarded and restored, `legend(ax)` not `legend()`,
   `mustBeMember` on every string option).
   *Verify:* draws into a caller-supplied axes without disturbing gca or the caller's
   hold state.
3. `analysis/plotSupernatantHeight.m` — stitches the legs of each `fld2x_*` chain the
   way `recreateFigsManriquez2026.m/loadChain` does (do not diverge) and writes
   `rec_supernatant_height.{pdf,png}` into `analysis/results/figures/recreation/`
   via the same `saveBoth` convention. Horizon and chain tags stated in the figure.
   *Verify:* the figure's final values equal the Python numbers in the progress file.
4. Report the numbers and whether winter or summer is taller.

## Files touched

- `src/@Results/getSupernatantHeight.m` (new)
- `src/@Results/plotSupernatantHeight.m` (new)
- `analysis/plotSupernatantHeight.m` (new)
- `analysis/results/figures/recreation/rec_supernatant_height.{pdf,png}` (new)

Nothing existing is modified. No simulation is started; no cosmos submission.
