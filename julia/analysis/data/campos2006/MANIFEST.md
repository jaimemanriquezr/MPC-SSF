# Campos, Smith & Graham (2006) — Extracted Data Manifest

Source paper: Campos, L. C., Smith, S. R., and Graham, N. J. D. (2006).
"Deterministic-Based Model of Slow Sand Filtration. I: Model Development."
*Journal of Environmental Engineering*, 132(8), 872–886.
File: `/Users/jaime/Research/SSF/documents/Campos2006.pdf` (16 pp).

All CSVs were transcribed directly from the PDF (pdftotext -layout + high-res
image rendering for verification). Every file has a header row, a `source`
column identifying the origin table/page, and a `units` column (or units in the
header) where applicable. Scientific notation normalized to E-notation
(e.g. `6.65E-04`). Numbers are transcribed exactly as printed — see Flags below.

## CSV files

### table1_filtration_runs.csv  (Table 1, p. 878)
Characteristics of the pilot-scale filtration runs used for calibration/verification.
Columns: source, filter (SSF A / SSF B), run, used_for (Calibration/Verification),
date_of_operation, duration_days, duration_hours, bed_depth_m.
5 runs (SSF A Run 1; SSF B Runs 1–4).

### table2_physical_parameters.csv  (Table 2, p. 880)
Physical parameters, stoichiometric coefficients, and constants (fixed inputs).
Columns: source, description, symbol, value, units.
Includes the **13 temperature correction factors θ** (theta_kga=1.066,
theta_kgb=theta_kgp=1.047, all others=1.08). Note: θ factors were assumed based
on Chapra (1997).

### table3_initial_parameter_ranges.csv  (Table 3, pp. 881)
Initial (prior) ranges of the 36 SSF model parameters for Monte Carlo calibration.
Columns: source, description, parameter, units, minimum, maximum, average, reference.
This is the primary source of rate constants and half-saturation constants with
units and ranges. (Table spans a column break across p. 881; both halves merged.)

### table4_best_ten_calibrated.csv  (Table 4, p. 882)
Best ten calibrated parameter sets (RB1–RB10) + Average, for Run 1 and Run 2 of SSF B.
Wide format: one row per parameter set, Run 1/Run 2 columns for each of
E0s, Ain, Eh, lambdaS0, ksa, krp, alpha. Units are given in Table 3.
Note: header uses **Eh** here but the same quantity is called **Fh** in Table 5.

### table5_calibration_before_after.csv  (Table 5, p. 882)
Six headloss parameters: initial range vs. calibrated range (min–max of the ten
best) for Run 1 and Run 2 of SSF B.
Columns: source, parameter, initial_range_min/max, source_reference,
calibrated_run1_min/max, calibrated_run2_min/max.

### table6_verification_strategies.csv  (Table 6, p. 883)
Three temperature-correction strategies for verifying Runs 1/2 in subsequent runs.
Columns: source, strategy_number, characteristic.

### intext_observations.csv  (running text, pp. 7, 9, 12)
Numeric experimental/derived values reported in prose (not in any table):
schmutzdecke thickness & biomass (SB) averages for Runs 1 & 2, adopted influent
protozoa concentration range, kgmaxb reduction in Strategy 3, SSF A Run 1
duration, and the numerical time step. Columns: source, description, run, value,
units, notes.

## Flags / interpretation notes

- **table4 RB2 krp_run1 = 7.90E-03**: transcribed exactly as printed (verified
  against a 300-dpi image crop). This value lies OUTSIDE the Run-1 calibrated
  range reported in Table 5 (6.6E-04 – 2.5E-03) and is ~4× the column average
  (1.95E-03). Likely a typographical error in the original paper (probably meant
  7.90E-04). Kept as printed; flagged here.
- **Eh vs Fh**: Table 4 labels the grazing-efficiency-related column "Eh"; Table 5
  labels the same parameter "Fh" (range 0.4–0.8). Both preserved as printed.
- Table 2 symbol `theta_Cga` is the temperature factor "for Cgz" per the
  description (apparent minor symbol/subscript inconsistency in the source).
- `ksCd` (Table 3) is written "kscd"/"ksCd" inconsistently in the paper; used ksCd.
- Some Table 3 parameters (a1, a2, c1, c2, kam, kew, ksl, xlm, fd, En, Yb, Yp) are
  empirical/light/grazing constants not all used downstream — included for completeness.

## Figure-only data (NOT digitized — no exact values in paper)

The following contain experimental data only as plotted curves; exact numbers are
not tabulated and were NOT digitized (per instructions, not fabricated):

- **Fig. 7 (p. 880)**: Sensitivity coefficients for uncovered SSF B (bar/plot of
  the ten most important parameters). Ranking is discussed but numeric SC values
  not tabulated.
- **Fig. 8 (p. 882)**: Four best calibration simulations of headloss vs time,
  Run 1 (a) and Run 2 (b), SSF B. Headloss-vs-time curves — figure only.
- **Fig. 9 (p. 882)**: Observed vs calculated headloss, Run 1 (a) / Run 2 (b) SSF B.
- **Fig. 10 (p. 883)**: Predicted headloss in Run 3, uncovered SSF B (Strategy 3).
- **Fig. 11 (p. 883)**: Verification of parameter set RB1 (Run 1) applied to SSF A.
- **Figs. 1–6**: Conceptual/schematic diagrams (schmutzdecke structure, filter
  representation, algae dynamics, microbial interactions, program flowchart) —
  no experimental data.

## Data NOT present in this paper

- No tabulated influent/effluent water-quality time series (turbidity, DOC, N, P,
  chlorophyll) — referenced as "extensive data sets" from Thames Water but not
  reproduced as tables here.
- No tabulated headloss-vs-time or biomass-vs-depth profiles (figure-only, above).
- No sand grain-size distribution table (only bed depth in Table 1 and porosity
  ranges E0/E0s in Table 3).
