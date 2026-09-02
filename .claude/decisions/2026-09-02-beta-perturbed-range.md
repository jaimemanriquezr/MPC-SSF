# β reported as PERTURBED in the sensitivity section; no re-run needed (2026-09-02)

Context: closes the sensitivity-section TODO on `beta_porosity` (β). The range was verified against the campaign
source and the raw measures by Claude on 2026-09-02, not by an independent check; Jaime to
confirm the wording before it feeds the manuscript.

## The decision

β (`beta_porosity`) is reported as **PERTURBED** in the manuscript sensitivity section, the
TODO on it is removed, and **no re-run is needed** — the OAT campaign already respected the
β ≥ 0.8 floor.

## Why

- `analysis/logOatCampaign.m:382-386` implements the 2026-08-18 constrained scheme:

  ```matlab
  beta = P("beta_porosity", "biofilm", 0.0, "Lund beta=0.99; constrained to [0.95, 0.99]", ...
           {@(f,m,gap) setModelArg(f, m, "BiofilmPorosity", 1 - gap)});
  beta.plus = 0.02;  beta.minus = 0.05;  beta.denom = log(0.05/0.02);
  ```

  It perturbs the **gap** (1 − β), not β itself. The two arms are β = 1 − 0.02 = 0.98 and
  β = 1 − 0.05 = 0.95 — both far above the 0.8 floor. There was never a ×2/×½ multiplier
  applied to β directly, unlike every other parameter in the campaign.
- The secant denominator is `log(0.05/0.02)` = ln 2.5, not the 2·ln 2 used to form the
  log-sensitivity for every other parameter in `logOatCampaign.m` (the ×2/×½ arms give a
  denominator of ln 2 − ln(1/2) = 2 ln 2). β's log-sensitivity is computed on a different base
  from the rest of the tornado. This is deliberate — the constrained gap scheme was adopted
  precisely because a symmetric ×2/×½ on β itself would either leave the porosity range
  (β > 1) or cross into an uncited regime, per the 2026-08-18 decision it implements — and it
  should stay flagged so nobody reads the β bar against the same log-2 scale as the rest of the
  tornado.
- The three campaigns record the same arms consistently, all `flag = OK`, `clog_driver = 0`:

  | scenario | file | I_rms |
  |---|---|---|
  | pulse | `analysis/results/oat/log_oat_pulse_beta_porosity/measures.csv` | 0.03858493428451 |
  | startup | `analysis/results/oat/log_oat_startup_beta_porosity/measures.csv` | 0.00583737689743338 |
  | flowstep | `analysis/results/oat/log_oat_flowstep_beta_porosity/measures.csv` | 0.0277381412575603 |

  and the masked tables agree (`I_rms_masked`, same arms, resolution-band restricted):
  `analysis/results/oat/measures_masked.csv` → 0.0389783785167506 (pulse),
  `analysis/results/oat/measures_masked_startup.csv` → 0.00553726894718241,
  `analysis/results/oat/measures_masked_flowstep.csv` → 0.0271141339724729. Every one of these
  six rows reads `flag_plus=OK, flag_minus=OK, clog_driver=0` — no β arm, in any scenario,
  violated the 0.8 floor or drove the filter to clogging.
- `analysis/PARAMETERS.md:67` and `:154` document the scheme: β = 0.99 nominal, runs at
  β = 0.98 (gap 0.02) and β = 0.95 (gap 0.05), denominator ln 2.5, citing the 2026-08-18
  decision.

## Alternatives rejected, and why

- **Masking the β row in the manuscript table**: nothing to mask — no arm violated the floor,
  so there is no invalid data to exclude, unlike the resolution-based masking applied to
  parameters that move the front arrival (e.g. dispersivity).
- **Re-running β as a single cosmos task** (~30 min): unnecessary for the same reason — the
  existing arms (β = 0.98/0.95) are already inside the valid [0.8, 1) range the manuscript
  wants reported, so a re-run would reproduce numbers already on disk.

## Consequences

- The sensitivity-section TODO on β is removed; β is written up as PERTURBED alongside the
  other 28 parameters, with a footnote/caveat that its log-sensitivity denominator (ln 2.5) and
  perturbation direction (on the porosity gap, not β itself) differ from the rest of the
  tornado — see the "Why" note above.
- No campaign re-run, no change to `analysis/results/oat/` or the ranking/tornado artefacts
  already regenerated 2026-09-02.
