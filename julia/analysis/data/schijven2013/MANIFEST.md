# Schijven et al. (2013) — Extracted data manifest

**Paper:** Schijven, J.F., van den Berg, H.H.J.L., Colin, M., Dullemont, Y., Hijnen, W.A.M.,
Magic-Knezev, A., Oorthuizen, W.A., Wubbels, G. (2013). "A mathematical model for removal of
human pathogenic viruses and bacteria by slow sand filtration under variable operational
conditions." *Water Research* 47, 2592–2602.
Source PDF: `/Users/jaime/Research/SSF/documents/Schijven2013.pdf` (11 pages).

Organisms: bacteriophage **MS2** (model virus) and **E. coli WR1 = ECWR1** (model bacterium).
Locations/filters: W = Weesperkarspel, L = Leiduin (Waternet); D = DUNEA; G = Groningen.

All numeric values transcribed exactly from the paper. Rate coefficients are day^-1 unless noted.

---

## CSV files

### `table1_experimental_conditions.csv`
Source: **Table 1** (p. 2594). Operational conditions of the 13 filtration experiments.
Columns: `source, experiment, location, temperature_C, age_schmutzdecke_days,
filtration_rate_cm_per_h, grain_size_mm, filter_bed_depth_m, porosity, dispersivity_cm, notes`.
- Temperature in °C; age in days; filtration rate cm/h; grain size mm; depth m; porosity dimensionless.
- **Flag:** Table 1 header prints the dispersivity unit as "[cm-1]", but the transport equations
  (Sec 2.5) give dispersivity a_L in "[m-1]". The unit label is ambiguous/likely a typo; values
  transcribed verbatim (0.38–2.1). Treat unit as uncertain.
- Notes column carries seeding info (L1 MS2-only; L5 ECWR1-only; D-sand is iron-oxide rich/stickier).

### `table2_ms2_fitted_rate_coefficients.csv`
Source: **Table 2** (p. 2598). Two-site kinetic model parameters fitted to MS2 breakthrough
curves (Hydrus-1D), plus overall removal rate lambda and calculated log10 removal.
Columns: `source, organism, experiment, katt1_per_day, kdet1_per_day, katt2_per_day,
kdet2_per_day, mu_l_per_day, mu_s_per_day, lambda_per_day, log10_removal_Cx_C0`.
- katt/kdet = attachment/detachment rate coefficients (sites 1 and 2); mu_l, mu_s = inactivation
  (die-off) rate coefficients of free and attached microorganisms; lambda = overall removal rate
  coeff (eq 5). All [day^-1] except log10 removal (dimensionless log10 Cx/C0).
- Note: MS2 not seeded/absent for L5 (ECWR1-only), so no L5 row here (matches paper).

### `table3_ecwr1_fitted_rate_coefficients.csv`
Source: **Table 3** (p. 2598). Same columns/units as Table 2, for ECWR1.
- No L1 row (L1 was MS2-only), matching the paper.

### `table4_schmutzdecke_parameters.csv`
Source: **Table 4** (p. 2599). Calibrated logistic-growth Schmutzdecke parameters, common to
both organisms and all filters.
Columns: `source, parameter, mean, standard_deviation, units, description`.
- f0 = 2.1e-4 ± 4.3e-6 [m·°C] (scale factor; gains °C dimension in final model eq 10).
- f1 = 9.9e-2 ± 6.0e-3 [day^-1] (rate coefficient; enters model as alpha*f1).
- SDs are from leave-one-out cross-validation.

### `table4_sticking_efficiencies.csv`
Source: **Table 4** (p. 2599). Calibrated sticking (collision) efficiency alpha, organism- and
filter-specific. Dimensionless.
Columns: `source, organism, filter_group, sticking_efficiency_alpha_mean,
sticking_efficiency_alpha_sd, units, notes`.
- MS2: W+L 9.2e-3, D 4.4e-2, G 5.6e-3. ECWR1: W+L 1.5e-1, D 7.1e-1, G 5.2e-2.
- ECWR1 alpha is ~15–20x higher than MS2 on the same filter (bacteria stick more than MS2).
- D-sand highest alpha (iron-oxide rich). SDs from leave-one-out cross-validation.

### `text_summary_values.csv`
Source: text/abstract/Fig. 2 (pp. 2592, 2599). Numeric values stated in prose that define
sensitivity/validation bounds. Columns: `source, quantity, organism, value, units, page, notes`.
- Removal ranges: MS2 0.082–3.3 log10; ECWR1 0.94–4.5 log10 (abstract).
  **Flag:** abstract ECWR1 max = 4.5 (exp D1), but Table 3 lists a higher value 5.1 (exp D3).
- Cross-validation RMSE = 0.30 log; prediction accuracy ±0.6 log.
- Bland–Altman: mean diff (obs–pred) = −0.028, 2SD = 0.61.
- Flow-rate sensitivity (log removal increase from lowering filtration rate): MS2 +0.06 (W), +0.26 (L);
  ECWR1 +0.07 (W), +0.27 (L) — demonstrates removal is relatively insensitive to flow within tested range.
- MS2 properties: diameter 26 nm; isoelectric point 3.5–4.1.

### `model_equation_constants.csv`
Source: Eqs. 6, 7, 10 (pp. 2595, 2599). Fixed constants/structure of the SSF model equations,
for reference (not fitted data). Includes the temperature-dependent dynamic-viscosity relation
(eq 7, Voss & Provost 2010): mu = 2.414e-5 * 10^(247.8/(273+T-140)).

---

## Data NOT extracted (and why)

- **Fig. 1 (p. 2597) — breakthrough curves** (C vs day for all experiments, 8 panels). Figure-only,
  **not digitized**. Axes: y = C (pfp/ml for MS2, cfp/ml for ECWR1), log scale ~1e-3 to 1e6;
  x = Day, ranges 0–8/0–12 (W,L) and 0–5 (D,G). Underlying fitted parameters are in Tables 2–3.
- **Fig. 2 (p. 2599) — Bland–Altman plot.** Only summary stats digitized (mean −0.028, 2SD 0.61);
  individual points figure-only, not digitized.
- **Fig. 3 (p. 2600) — predicted MS2 removal vs temperature time series** (2006–2011, 358 temperature
  data points from full-scale Weesperkarspel/Leiduin). Figure-only, **not digitized**; underlying
  temperature dataset not tabulated in the paper.

## Notes / caveats
- Some Table 2/3 entries have katt2=kdet2=0 (experiments G1, G2) — one-site behaviour; transcribed as 0.
- mu_l = 0 and mu_s = 0.056 for MS2 W2; mu_l = 0.00 for ECWR1 W2 — transcribed exactly.
- lambda is a fitted/derived overall rate (eq 5), generally ≈ katt1 except low-temperature
  experiments W2, G2 where kdet1 is relatively high and lambda < katt1.
