# Reference methods digest: sensitivity analysis & calibration in prior SSF models

Prepared for the slow-sand-filtration (SSF) sensitivity-analysis writeup. Two source
papers are summarized below. Page citations refer to the PDF page numbers of each
article (Campos2006 = J. Environ. Eng. 132(8):872-886; Schijven2013 = Water Research
47:2592-2602). Quotes are kept short and marked; where a detail was not present or not
legible, this is stated explicitly.

---

## 1. Campos, Smith & Graham (2006) — Deterministic-Based Model of SSF, Part I

Mechanistic (deterministic) SSF model predicting headloss and filtration-coefficient
development, solved by finite-difference (Euler) with the sand bed discretized into 1-cm
sublayers (p. 878). Calibrated/verified against pilot-scale headloss data.

### Sensitivity analysis method
- Called an "A Priori Sensitivity Analysis": run *before* calibration to select the most
  sensitive parameters and so reduce the number of calibration simulations (p. 879).
- Scheme: strict **one-at-a-time (OAT)** — "holding all parameters constant while one is
  varied"; the paper explicitly notes this "ignores parameter interdependence" (p. 879).
- Model output examined: **headloss** response (p. 879). Data used: the first filtration
  run of uncovered filter SSF B (p. 879).
- Metric: a **sensitivity coefficient (SC)**, following Schütze et al. (2002), defined as
  the ratio of the coefficient of variation of the model output y to the coefficient of
  variation of the perturbed parameter p_x:
  - SC = COV_y / COV_px  (Eq. 31, p. 879).
  - COV_y = (100 / sqrt(n0-1)) * (s_y / mean(y)) — relative CV of output (Eq. 32).
  - COV_px defined analogously from the mean and std dev of the sampled parameter values
    (Eq. 34, p. 880). n0 = number of data points.
  - So SC is a dimensionless, variance-based normalized-sensitivity metric (see also the
    Notation list, p. 886: "SC = sensitivity coefficient (dimensionless)").
- Ranking: parameters ranked by SC magnitude; **Fig. 7** plots the SC for the ten most
  important parameters for uncovered SSF B (caption note: "only ten most important
  parameters are shown", p. 880). The individual bar values in Fig. 7 were not extracted
  as text (figure not OCR-legible); only the identity of the top parameters is given in
  the running text.

### Calibration method
- Calibration limited to the **six most sensitive** parameters identified above (p. 879).
  The six were: initial porosity of schmutzdecke (E0s); schmutzdecke growth rate (Ain);
  grazing efficiency of protozoa (Eh); initial filtration coefficient for schmutzdecke
  (lambda_S0); half-saturation constant for protozoa growth on algae (ksa); and protozoa
  respiration rate (krp) (pp. 879-880). All other parameters were fixed at literature /
  measured values (Tables 2-3, e.g., initial clean-bed porosity from personal
  communication; temperature-correction factors assumed from Chapra 1997, p. 880).
- Objective function: **least-squares / sum of squared errors**, unweighted (all measured
  values treated as equally important, no weighting):
  - SSE = sum_i (y_obs,i − y_i)^2  (Eq. 35, p. 880).
- Optimization approach: **Monte Carlo** search — 300 Monte Carlo simulations using
  **uniform random sampling** across specified initial parameter ranges; the ten best
  parameter sets (lowest SSE) were retained per run (p. 880). Initial ranges vs. calibrated
  ranges compared in Table 5.
- Sequential-run handling: after calibrating filtration Run 1, the simulated absolute
  specific deposit remaining in the sand (after cleaning, 5 cm removed incl. schmutzdecke)
  was used as the input state to calibrate Run 2 (p. 879).
- Calibration/verification split: driven by pilot operating records (Table 1, p. 877).
  SSF B Run 1 = Calibration; SSF B Run 2 = Calibration; SSF B Runs 3 & 4 = Verification;
  SSF A Run 1 = Verification. Verification applied Run-1/Run-2 calibrated parameters to
  simulate Runs 3/4 and to the independent filter SSF A.

### Identifiability, correlation, seasonal/temperature dependence
- Identifiability/correlation: the OAT scheme "ignores parameter interdependence"
  (p. 879) — no formal correlation/identifiability analysis is reported. The authors note
  that only minor differences between initial and calibrated ranges indicated that
  "several combinations of parameters gave good predictions of headloss" (p. 880),
  i.e., an implicit acknowledgement of non-uniqueness / equifinality.
- Seasonal/temperature dependence: emphasized. Calibrated parameters (notably schmutzdecke
  growth rate Ain, and headloss empirical parameters c1, c2, b) were found to vary between
  runs, attributed to seasonal/environmental factors (Run 2 in winter → smaller Ain range)
  (pp. 880, 883). Verification used an explicit **temperature-correction** of headloss
  parameters: P_v = (P_c / TM_c) * TM_v, scaling calibrated parameters by the ratio of
  average run temperatures (Eq. 36, p. 883). Model verification "suggested [parameters]
  may be influenced by seasonal effects" (p. 885).

### Experimental design (data source)
- Pilot plant: Kempton Park Advanced Water Treatment (AWT), Thames Water Utilities, treating
  lowland reservoir water; three parallel slow sand filters (two uncovered, one covered);
  modeling used the two uncovered filters SSF A and SSF B (pp. 877-878).
- Runs: four filtration runs of SSF B plus first run of SSF A; run durations 99-152 days;
  bed depths 0.59-0.70 m (Table 1, p. 877).
- Measured quantities: primarily **headloss and flow rate** over consecutive runs; bulk
  specific deposit was NOT measured (its range bounded by initial medium porosity), and
  protozoa counts were not routinely available (pp. 878). Calibration therefore targeted
  headloss.

---

## 2. Schijven et al. (2013) — Mathematical model for virus/bacteria removal by SSF

Fits a **two-site kinetic attachment/detachment model** (Hydrus-1D) to pilot-column
breakthrough curves, then builds a steady-state log-removal predictor that adds a
Schmutzdecke term, and studies removal vs. operational conditions.

### Sensitivity analysis method
- Purpose (Section 2.10, p. 2596): demonstrate the model under variable operational
  conditions and "investigate the sensitivity of the predicted removal" of MS2 and E. coli
  WR1 to operational conditions.
- Scheme: **scenario-based forcing sweeps** rather than a formal parameter-perturbation /
  SC-coefficient analysis. Removal of MS2 and ECWR1 was predicted for two full-scale
  locations (Weesperkarspel and Leiduin) driven by **358 full-scale water-temperature
  records (2006-2011)**; a scenario with a one-year-old Schmutzdecke scraped off at low
  temperature (2007) and high temperature (2009); and each location run at **two filtration
  rates** (pp. 2596, 2599). The differing grain size between the two locations acts as an
  additional varied condition.
- Metric: no single scalar sensitivity coefficient is defined. Sensitivity is expressed as
  the change in predicted **log10 removal** when a condition changes. Reported effects
  (p. 2599): lowering filtration rate from 45→30 cm/h (Weesperkarspel) and 30→22 cm/h
  (Leiduin) raised MS2 removal by only 0.06 and 0.26 log10 (ECWR1: 0.07, 0.27), whereas
  temperature and Schmutzdecke age drove large periodic variability. Scraping the
  Schmutzdecke caused the largest removal drop when done at high temperature.
- Ranking (qualitative): predicted removal is **strongly sensitive to water temperature
  and Schmutzdecke age**, and only **weakly sensitive to filtration rate** within the tested
  range (pp. 2599, 2600 Abstract/Discussion). Grain size explains the between-location
  removal difference.

### Calibration method
Two nested fitting stages:
1. **Breakthrough-curve fitting (transport parameters).** All salt-tracer and microorganism
   breakthrough curves were fitted with **Hydrus-1D**'s two-site kinetic model (Simunek et
   al. 2005) (p. 2595). Salt tracer → porosity n and dispersivity/dispersivity alpha_L
   (Table 1). Microorganism curves → the two-site rate coefficients per experiment:
   k_att1, k_det1, k_att2, k_det2, plus inactivation rates mu_l, mu_s (Tables 2-3, p. 2597).
   The overall removal-rate coefficient lambda is derived (Eq. 5). Inactivation rates
   mu were measured **independently** (batch die-off in supernatant, linear regression of
   log-transformed counts, Section 2.4, p. 2595), i.e., fixed inputs, not free fit params.
2. **Schmutzdecke-model calibration (removal-predictor parameters).** A logistic growth
   term f0*(1 − exp(−f1*a)) with age a was added on top of colloid-filtration attachment
   (sticking efficiency alpha) (Eqs. 6-8, pp. 2595-2596). The three parameters **alpha, f0,
   f1** were estimated by fitting the log-removal objective function
   ln(C/C0)[a_i, f0,i, f1,i, dp,i, dc,i, n_i, u_i, x_i] to each experiment's removal data.
   - Objective function: **maximum likelihood** — minimize a log-likelihood function
     (McCullagh & Nelder 1989; Eq. 9, p. 2596) by numerical optimization in
     **Mathematica 8.0.4** (Wolfram).
   - Parameters fixed vs. fitted: operational conditions (dp = organism size, dc = grain
     size, u = filtration rate, T = temperature, n = porosity) are fixed covariates; alpha,
     f0, f1 are fitted. alpha may differ per organism/filter; f0, f1 may be shared across
     filters/organisms — tested by **likelihood-ratio tests** (Section 2.8, p. 2596): the
     model with fewer distinct a/f0/f1 values is preferred if the LR is not significant.
- Calibration/validation split: **leave-one-out cross-validation** (24 log-removal data
  points from MS2 + ECWR1 experiments). Refit on 23, predict the excluded experiment,
  repeated 24 times; report **RMSE** of predicted vs. observed log removal, plus a
  **Bland-Altman** limits-of-agreement plot (Section 2.9, p. 2596; RMSE = 0.30 and
  agreement within ~0.6 log, Table 4 / p. 2599).

### Identifiability, correlation, seasonal/temperature dependence
- Identifiability handled via **likelihood-ratio parsimony tests** to decide whether a, f0,
  f1 can be shared across filters/organisms (fewer parameters preferred when LR not
  significant) (p. 2596). f0, f1 were found common (single mean values, Table 4); alpha
  differed by filter/location and organism (Table 4, p. 2599).
- Parameter correlation among the two-site rates is discussed qualitatively: lambda ≈ k_att1
  because k_det1, mu_l, mu_s are small and k_det2 large (Eq. 5; pp. 2595, 2597) — i.e.,
  attachment to type-1 sites dominates removal, so the other rates are weakly constrained.
- Temperature dependence: central. Removal depends strongly on temperature (via dynamic
  viscosity and diffusion in colloid-filtration theory, and via Schmutzdecke activity);
  the sensitivity study is driven by multi-year temperature records (pp. 2597, 2599-2600).
  The Schmutzdecke term makes removal age-dependent, coupling to seasonal scraping timing.

### Experimental design (data source)
- Pilot-scale slow-sand filters at three Dutch utilities: Weesperkarspel (W), Leiduin
  (L, incl. DUNEA sand filter D) and Groningen (G). Filter surfaces ~1.6 m × 1.6 m
  (Leiduin) and a 1.25-m-diameter, 1.25-m-deep pilot (Groningen); filter-bed depths ~1-1.5
  m with ~1.2 m supernatant (pp. 2593-2594). 24 experiments (W1-W3, L1-L5, D1-D3, G1-G2,
  etc.; Table 1).
- Tracers/organisms: NaCl salt tracer (550 mg/l, 24 h seeding) for pore-velocity and
  dispersivity; bacteriophage **MS2** (26 nm, conservative model virus) and **E. coli WR1**
  (ECWR1) seeded into supernatant, ~30 effluent samples over 5-7 days per run
  (Sections 2.3-2.4, pp. 2594-2595).
- Varied operational conditions (Table 1, p. 2594): water temperature, Schmutzdecke age
  [days], filtration rate [cm/h], grain size [mm], filter-bed depth [m], porosity,
  dispersivity. Measured output: log10 removal (Cx/C0) of MS2 and ECWR1.

---

## Notes on legibility
- Campos2006 Fig. 7 (the SC bar chart) is a figure image and its per-parameter SC values
  were not OCR-recoverable from `pdftotext`; only the identity/ordering of the top
  parameters given in the running text (pp. 879-880) was used. If exact SC values are
  needed for the writeup, render p. 880 as an image.
- Campos Tables 2-5 (parameter ranges/values) came through with column layout scrambled;
  individual numeric ranges were not transcribed here (not needed for the methodology
  digest) but are legible enough to extract on demand from pp. 879-882.
- Schijven2013 equations 1-9 rendered as garbled math in text extraction but the surrounding
  prose defining each method (objective functions, LR tests, cross-validation, sensitivity
  scenarios) was fully legible; no methodological detail was invented.
