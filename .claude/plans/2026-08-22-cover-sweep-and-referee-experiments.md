# Experiment plan: answering the Cleaner Water referees, and making covered/uncovered bite

## Context

`Manriquez2026` (CLWAT-D-26-00074) came back **major revision**. Two referees; the
experiment-bearing comments are R1.2 (no verification against published data), R1.4 (the
oxygen anomaly), R1.6 (what mechanism actually removes the pathogen), R2.2 (which
parameters matter) and R2.3 (soften "predicting" to "qualitatively reproduces"). R1.3,
R1.5, R2.1 and R2.4–2.6 need text, not runs, and are **out of scope here** (scope decision:
experiments only).

Two things set the shape of this plan.

**First**, the literature synthesis (`reports/literature-synthesis-2026-08.typ` §6) recast the
model/field comparison in terms of *responses to perturbations* rather than levels. Responses
are ratios between paired runs, so they divide out β, ρ_P, the attachment prefactor and the
absolute influent scale — all of which are unresolved and all of which block level
comparisons. Every experiment below is therefore specified as a **paired contrast reported as
a ratio**, which additionally makes it robust to the open mass-conservation defect
(`.claude/decisions/2026-08-20-fixed-influent-flip-evaluation.md`), since both arms share the
transport operator and a defect there is common-mode.

**Second**, the covered/uncovered experiment (E2) produces a ~4 % biomass difference where
Campos2002 measured 4× at full scale, and under the corrected stack it is *worse*
(covered 0.342 vs uncovered 0.351, ratio 0.97 — `reports/results-redo-poc-2026-08.typ`). That
is a qualitative failure, not a magnitude one, and it needs diagnosis before it can be
reported. Exploration turned up **four** candidate blockers, not one, and they are separable.

**Model of record for everything below: the corrected stack** (audit kinetics
μ_HET = 1.008 /d, μ_PHO = 3.0 /d, source-faithful half-sats, endogenous phototroph
respiration k_ra = 0.276 /d, field-strength influent). Manuscript-stack cells are retained
only as continuity controls where noted.

---

## Part 0 — What is already answered and needs no runs

Check this before building anything.

| Comment | Status | Evidence |
|---|---|---|
| **R2.2** — which parameters are most influential, naming attachment/detachment, biofilm excess velocity, light attenuation, PAT inactivation/bacterivory | **Already done.** All four groups are in the completed log-OAT campaign: `attach_sand`, `sand_pathogen`, `detach_scale`, `transport_P`; `kappa`, `zeta_0`, `zeta_1`; `light_att_water`, `light_att_sand`, `attenuation_P`; `bacterivory`, `inactivation`, `marker_growth`. | `reports/log-oat-report-2026-08.typ` (26 params × 3 scenarios, 55 runs each, all `flag = OK`); `reports/sensitivity-campaign-2026-08.typ` (three tiers complete) |
| **R1.4** — explain the oxygen increase | **Mechanism identified and reversed.** Corrected model gives −3.7 mg/L against a 9.10 influent, bed minimum 3.34. | `reports/kinetics-audit-2026-08.typ`; `reports/results-redo-poc-2026-08.typ` |

Two notes that must go into the rebuttal rather than into new runs:

- The campaign finds **light parameters negligible in every scenario**, which contradicts
  R2.2's stated expectation that light attenuation coefficients "are likely to have
  substantial effects". Say so directly — it is a stronger answer than agreeing.
- The OAT campaign matured for **3 days at 30 cells**, not 30 days at 100. A referee may
  object that the sensitivities were measured on an immature filter. Decide whether to
  disclose-and-defend or to re-run one scenario at manuscript resolution (~55 runs × 8 min
  ≈ 7 core-hours; affordable but not free).

---

## Part 1 — New experiments

Four experiments. All are paired contrasts. Durations assume 100 cells, `AdaptiveMaxDt=3e-6`,
measured at **≈ 0.8 min per simulated day** (30 d ≈ 25 min, 10 d ≈ 8 min, 104 d ≈ 85 min).

### X-BIO — Bellamy's biological-activity gradient → R1.6, R1.2

The only experiment in the corpus that isolates the biological contribution with the physics
fixed. Three arms on a 30-day mature filter under constant marker feed, varying only the
biological removal channels:

| Arm | Change | Bellamy target |
|---|---|---|
| `suppressed` | `BacterivoryRate` → 0, `InactivationRate` → 0 | 0.4 log (chlorinated filter, 60.1 %) |
| `reference` | corrected-stack nominals | 1.6 log (control, 97.5 %) |
| `augmented` | influent DOM raised until the bed draws 3–4 mg/L O₂ across it, matching Bellamy's synthetic-sewage loading | 3.0 log (nutrient-fed, 99.9 %) |

**Acceptance:** monotone increase, ≥ 1.5 log spread between `suppressed` and `augmented`.
**Cost:** 3 × 10 d = 24 min + shared mature cache.
**Why first:** answers R1.6 with a number rather than a paragraph, needs no new physics, and
the `augmented` arm requires a short calibration loop on DOM (2–3 throwaway 2-day runs) to hit
the 3–4 mg/L drawdown — budget for that.

### X-TSPLIT — the biological/physical temperature split → R1.6, R1.2, R2.3

The single structural claim every experimental paper in the corpus agrees on (synthesis §4.3):
cold suppresses biological removal and leaves physical removal untouched. A 2×2,
{19 °C, 3 °C} × {biology on, biology off}, reported as **two temperature differences**.

**Targets:** Δ_T(biology on) ≈ 0.8–2.1 log (the seeded-indicator band, 0.05–0.13 log/°C over
16 °C, corrected synthesis §4.2); Δ_T(biology off) ≈ 0.
**Acceptance:** Δ_T(on) ≥ 5 × Δ_T(off), physical arm flat to within 0.1 log.
**Cost:** 4 × 10 d = 32 min, **plus a winter mature cache** (30 d at 3 °C ≈ 25 min) which does
not currently exist.
**Blocker:** `runMarker` hardcodes `theFilter(opts, 19, lightSummer())` at
`analysis/manuscriptExperiments.m:337`, and `getMature` only ever builds a summer cache.
**Note:** if the physical arm shows appreciable temperature sensitivity, the cause is an
Arrhenius factor on a rate that should not carry one — a finding in itself, worth reporting.

### X-INFL — influent-strength response → R1.2, R2.2

Sweep influent PHO from field strength to Table B.1 and beyond (5e-4 → 1e-2 → 1e-1 kg/m³),
everything else fixed, reporting φ_b, its depth distribution, effluent O₂ and marker log
removal as *functions* of influent strength.

**Target:** Bae2023 — a 39 % higher particle load in summer changed removal not at all
(99.3 vs 98.7 %). Removal should be insensitive; standing biomass should not be.
**Acceptance:** log removal flat to within 0.2 log; φ_b monotonically increasing.
**Cost:** subsumed into X-COVER Stage B (the `I` factor). Do not run separately.
**Why:** the audit's central boundary-condition recommendation currently rests on one probe at
one value; a response curve turns it into a characterised behaviour.

### X-COVER — the covered/uncovered driver sweep → R1.2, R2.3

Full design in Part 2.

### Not in scope, but flag in the limitations

Grain-size dependence (Bellamy 0.82 log over 0.128–0.615 mm; Schijven ∝ 1/d_c) **cannot be
tested** — there is no `d_c`, collector efficiency or Tufenkji correlation anywhere in `src/`.
Same for inert particles, hydraulic feedback and EPS. These are structural absences, not
calibration failures; no run will produce them.

---

## Part 2 — The covered/uncovered driver sweep

### The four candidate blockers

Campos2002 measured, on two full-scale Thames beds ripened in parallel from bare sand:
uncovered 0–2 cm 60 µgC/g (126.8 peak) against covered 15; 0–10 cm 22.5 against 5.6;
**no schmutzdecke at all** on the covered bed; and **identical TOC/DOC removal** (25/23 %
vs 23/23 %). Our model gives a 4 % difference, 0.97 under the corrected stack.

1. **Light amplitude.** The cover is a scalar on `f.LightIrradiation` (`e2_covered`,
   `manuscriptExperiments.m:287`), applied before `exp(-eta)`, so a factor s is exactly
   `-ln s` added to η uniformly — amplitude only, no shape change. Code uses 1 %; the
   manuscript text says 0.1 % (a documented deliberate divergence, `analysis/PARAMETERS.md:124`).
2. **Nutrient limitation.** Influent HPO₄ = 0 exactly, so phototroph growth runs on P recycled
   from death and respiration — biomass-proportional and light-independent. A 4× light cut
   moved winter growth 6 % (`.claude/decisions/2026-08-20-winter-summer-flip.md`).
   *Complication:* influent NH₄ = 2.0e-5 against the published `K_NH4 = 1.2e-2` is a Monod term
   of 0.0017, a tighter throttle than P. Under the corrected stack (`K_NH4 = 2.0e-5`) NH₄ sits
   near 0.5 and P binds instead. **The limiting nutrient differs between the two kinetic
   stacks**, so a P-only sweep would answer the wrong question — hence the 2×2 add-back below.
3. **Influent dominance.** Table B.1 PHO = 1e-2 is 15–20× field. At 7.2 m/d that delivers
   2.16 kg/m² over 30 d against a standing stock O(0.1) kg/m² — trapping outruns in-filter
   growth ~20:1, so covering cannot matter regardless of light or nutrients.
4. **Euphotic depth (found during design; not previously on the list).**
   `LightAttenuationCoeffSand = 1500` applied to `(1−ε₀)(δ/2 + z)` is an effective 900 m⁻¹ in
   the bed. At `addGridPoints(100)`, dz ≈ 9.95 mm, so the first submerged bed cell already sits
   at e^−η ≈ 1.4e-5 and the second at ~1e-9. **The euphotic depth in sand is ~5 mm — half of one
   grid cell.** Essentially all photosynthesis happens in the supernatant; bed biomass responds
   to light only second-hand, through settling and attachment. Campos's 0–2 cm layer is two
   cells, both effectively dark.

Blocker 4 also implies the reported "4 %" may be partly a **metric artefact**: φ_b integrated
over the domain is dominated by HET, POM and enclosed water, all light-independent. Stage A
resolves that for 25 minutes of wall time, and it must be resolved before anyone changes model
structure.

### Factors and levels

Shared and locked for every cell: T = 19 °C, `NCells=100`, `MaxDt=3e-6`, `Kappa=1e-7`,
`Zeta0=1e2`, `Zeta1=1e-2`, detachment `0.14*sqrt(|v|/18)`, `NormalizedLight=true`,
`ImplicitOsmosis=true`, `lightSummer`, **clean start** (matching Campos, whose beds were both
resanded and ripened in parallel), hourly frames, corrected kinetic stack.

| id | factor | levels |
|---|---|---|
| **S** | cover multiplier | 1, 0.1, 0.01, 0.001, 0 — uniformly spaced in η (0, 2.3, 4.6, 6.9, ∞). Includes the code's 1 %, the text's 0.1 %, and 0 as the physical bound of a light-excluding membrane |
| **N** | influent nutrient add-back | none / +P (HPO₄ = 5.0e-4) / +N (NH₄ = 1.0e-2) / +NP |
| **I** | influent PHO | 1.0e-2 (Table B.1) / 5.0e-4 (Campos2006b field) |
| **K** | kinetic stack | corrected throughout; one manuscript-stack pair retained as a continuity control |

### Staging, with gates

Gates are **human** — one look at `pairs.csv` between stages. Automating them would spend
Stages B and C blindly, which is what staging exists to prevent.

**Stage A — light ladder + metric audit.** 5 runs × 30 d ≈ 25 min wall.
S ∈ {1, 0.1, 0.01, 0.001, 0}, all else nominal. Yields the elasticity ε_L = d ln Q / d ln s per
QoI and the S = 0 ceiling.
- *G-A1* — if `R_B_0_2` at S = 0 is ≥ 4× while `R_phib_0_2` ≈ 1, the failure was a **metric
  artefact**. Fix the reported quantity, stop.
- *G-A2* — if `R_B_0_2` at S = 0 is < 1.5×, light amplitude is not the gate. Go to Stage B.
- *G-A3* — if ε_L collapses toward 0 at small s, something *compensates*; read the Part 3
  diagnostic to see which Monod term rises.

**Stage B — nutrient × influent, paired on cover.** 14 new runs × 30 d (16-cell grid minus the
two already run in Stage A).
S ∈ {1, 0.001} × N ∈ {none,+P,+N,+NP} × I ∈ {1e-2, 5e-4}.
- *G-B1* — moves only when N ≠ none → **nutrient limitation gates**; the +P/+N/+NP contrast names
  which, and the argmin diagnostic confirms independently.
- *G-B2* — moves only at I = 5e-4 → **influent dominance gates**.
- *G-B3* — moves only at both → Table B.1's influent is the single change that unlocks E2.
- *G-B4* — nothing moves → **euphotic depth gates**. Trigger Stage B′.

**Stage B′ — conditional on G-B4 only.** 4 runs × 10 d. Two counterfactual pairs, neither a
proposed model change: `LightAttenuationCoeffSand` 1500 → 150 (euphotic depth 5 mm → 5 cm), and
`NCells` 100 → 300 (dz ≈ 3.3 mm, resolving the euphotic layer). These convert "we think it is
the euphotic depth" into a demonstration. Budget the N = 300 arm generously; it may need a
reduced `MaxDt`.

**Stage C — confirmation at Campos's own duration.** 4 runs × 104 d ≈ 85 min.
Covered/uncovered at the winning configuration, plus covered/uncovered at baseline as the null.
104 d rather than 90 removes a duration confound for 14 min per run.

### Quantities of interest

Campos reports **µg C per g dry sand**. Conversion from our state:
ρ_sand = (1−ε₀)·2650 = 1590 kg/m³; COD→C for C₅H₇O₂N = 0.531/1.42 = 0.374 kgC/kgCOD; hence

> **κ = 235.2 (µgC/g dry sand) per (kg COD/m³)**, and
> `B̄(0→d) = κ · (1/d) · ∫₀^d [X_HET(z) + X_PHO(z)] dz` for d = 0.02, 0.10 m

Sanity check: Campos's 60 µgC/g back-maps to X_att = 0.255 kg COD/m³ — right order for a
ripened top layer. Because Campos normalises *per gram of sand* and ε is constant inside the
bed, the correct depth average is the **plain arithmetic mean**, not the ε-weighted one; the
ε-weighting rule applies to conserved-mass integrals, which these are not. Report both and
label which is which.

Per cell record: the Campos targets (`B_0_2`, `B_0_10`, ± POM variants), the PHO-only
components, peak and its depth, the E2-continuity `phib_*`, a **schmutzdecke presence/absence
QoI** (`mat_mass` over z<0 and the extent of {z<0 : φ_b>0.05} — "no schmutzdecke at all" is a
presence measurement and needs its own metric), **DOC/TOC removal** (Campos's equally
load-bearing second result), the limitation reductions from Part 3, the photosynthesis
partition across supernatant / 0–2 cm / 2–100 cm, and `flag`/`tFinal`/mass residuals.

Reported result is the **covered ÷ uncovered ratio** per matched configuration:
`R_B_0_2` (target 0.125–0.25), `R_B_0_10` (0.25), `R_mat` (≈ 0), `ΔDOCrem`/`ΔTOCrem` in
percentage points (|Δ| < 3), and `eps_L = ln R / ln s`.

**Joint success criterion** — `R_B_0_2 ∈ [0.12, 0.25]` **AND** `|ΔTOCrem| < 3 pp` **AND**
`R_mat < 0.1`. A configuration that wins the biomass ratio by also destroying removal has not
reproduced Campos.

**Mass-defect guard:** record the ε-weighted mass-balance residual per arm and **reject any
pair whose two residuals differ by more than 5 %**. Common-mode cancellation is the argument
for trusting ratios; make it auditable rather than assumed.

---

## Part 3 — Instrument the Liebig min (do this first)

The limiting nutrient should be **measured, not inferred from sweep outcomes**. The Liebig min
is computed in exactly one place and its argmin is already available and thrown away.

**`src/@State/simulate.m:608-621`, `evaluateReactions`** — return the argmin alongside:

```matlab
function [rx, monod, lim] = evaluateReactions(local, K, phi, mu, I, orders)
    monod = ones(size(phi,1), size(mu,2));
    lim   = zeros(size(phi,1), size(mu,2), "uint8");   % 0 = no Monod terms
    for i = 1:size(monod,2)
        ii = K{1,i};
        if any(ii)
            k = K{2,i};
            mon_term = (local(:,ii) + realmin)./(k + local(:,ii) + realmin);
            [monod(:,i), j] = min(mon_term, [], 2);
            gidx = uint8(find(ii));
            lim(:,i) = gidx(j);
        end
    end
    rx = phi.*mu.*I.*monod.*product;   % bit-identical
end
```

Three additive edit sites, all behind a new `parameters.RecordLimitation (1,1) logical = false`
in the `arguments` block so every golden stays byte-identical:

1. `:608-621` as above (`gidx` indexes `[particles, liquids, quotients]`; build the name map
   once at `:49-63` where `quotientNumIdx` is already in scope).
2. `:346` and `:348` — capture the extra outputs for the **biofilm and flowing** regions. Both
   matter: photosynthesis is in the supernatant flowing phase, standing stock is in the biofilm
   phase. `:347` (enclosed) can stay.
3. Frame recording — preallocate beside `concFramesBiofilm` (`:150-153`), write in the frame
   block (`:557-566`), export at `:597-602` as `results.Frames.Limitation`. Store `lim` as
   `uint8`, `monod` as `single`. Also record `lightAttenuated` and `lightFactor(:,2)`, both
   already in scope at `:322-330`.

*Wrinkle:* the pre-loop frame write at `:165-171` runs before `monodB` exists. Initialise to
NaN/0 and let frame 1 carry them; do not restructure the loop.

Cost: no extra floating-point work, ~13 MB per run, zero when the flag is off.

**What it buys.** Growth = φ·μ(T)·**l**·**min_j m_j**·X, and a cover scales `l` by s *linearly*.
Any sub-linear response therefore lives entirely in the feedback, which is visible only in
min_j m_j. So the decisive statistic is the elasticity decomposition

> Δln(growth)/Δln(s) = Δln(l)/Δln(s) + Δln(m\*)/Δln(s)

The first term is ≈ 1 by construction (with a measurable Steele-concavity correction). If the
total is ≈ 0 then Δln(m\*)/Δln(s) ≈ −1: the limiting substrate recovers exactly as fast as light
is removed. **That is the falsifiable signature of the recycling loop, measured rather than
assumed.** Recording the argmin identity additionally settles whether the manuscript stack is
N-limited while the corrected stack is P-limited — which decides whether Stage B's add-back
levels relieve the right nutrient, and is why this lands *before* Stage B.

---

## Part 4 — Infrastructure to build

No generic sweep driver and no results collector exist. Follow `phoBloomStudy.m` (struct-array
of named configurations) and `probeSeason90.m` (the modern probe template: relative paths via
`fileparts(mfilename("fullpath"))`, an `arguments` block, ε-weighted integrals via
`computePorosity`). **Do not** copy the 16 legacy probes — they hardcode a stale worktree path
missing the `1d/` segment and cannot run as written.

| File | Purpose |
|---|---|
| `analysis/probes/coverCases.m` | **New.** `cases = coverCases(stage)` — the registry and single source of truth. Fields: `stage, tag, cover, phoIn, nh4In, hpo4In, tsim, ncells, maxdt, kinetics, etaSand, pairKey`. `pairKey` is the tag minus the cover field, so the collector pairs on it and never parses filenames. One `coverTag(cfg)` helper used everywhere. |
| `analysis/probes/probeCover.m` | **New.** `probeCover(idx, opts)` — the cell runner. Builds filter + model from `cases(idx)`, `simulate(..., RecordLimitation=true)`, computes the Part 2 QoIs, saves `data/cover_<tag>.mat` and `results/cover_sweep/profile_<tag>.csv`. |
| `analysis/collectCover.m` | **New.** `[T,P] = collectCover(stage)` — the collector. **Assert `isequal(loaded.cfg, cases(i))`** so a stale `.mat` from an edited registry is an error, never a silent reuse; this is the most important line in the file. Emits `summary.csv`, `pairs.csv`, `limitation.csv`. |
| `analysis/probes/probeBio.m` | **New.** X-BIO's three arms; reuses `runMarker`'s influent/log-removal idiom. |
| `analysis/probes/probeTSplit.m` | **New.** X-TSPLIT's 2×2. Needs a **winter mature cache**. |
| `slurm/cover_sweep.sbatch` | **New.** Indexed off the registry via `$SLURM_ARRAY_TASK_ID`, *not* parallel bash arrays — those are the maintenance trap in `season90.sbatch:12-15`. |
| `src/@State/simulate.m` | **Modify.** The three Part 3 edits, default-off. |

### The `getMature` cache-key trap

`manuscriptExperiments.m:190-201` keys on `NCells` + smoke flag only. It ignores `Respiration`,
`PGExcess`, `Kappa`, temperature, the light handle and the whole influent vector — so any run
that matures under varied conditions silently reuses a stale summer state.

Rules: **the cover sweep never calls `getMature`, `runSummer` or `manuscriptExperiments`** —
every cell is a clean start, which is both correct protocol and complete immunity. All output
paths are `coverTag(cfg)`-derived. The only sentinel is on the config-tagged path, always with
an `opts.Force` escape.

X-BIO and X-TSPLIT *do* need mature states, so they must either use a distinct `OutRoot` per
variant or extend the cache key. **File `getMature`'s key as its own bug** — it silently
corrupts E4/E5/E6/E10 across any kinetics change — and fix it independently rather than making
this plan depend on it.

### Execution environment

Per `probe-jobs-on-cosmos`: sweeps go to cosmos via `sbatch`, never local MATLAB batch.
Account `-A lu2026-2-100` (not `lu2025-7-124`, the HDG project's), partition `lu48`,
`module load matlab/2025b`, `#!/bin/bash -l` required, paths relative to the submit directory,
`slurm/logs/` must pre-exist. Deploy by `rsync` one-way; **the cosmos trees are stale copies
still carrying the retired `julia/` directory and the old wrong-account scripts, so re-rsync
before submitting anything**. Cosmos authentication requires Jaime's live ControlMaster socket
(`ssh -o ControlPath=~/.ssh/cm-cosmos cosmos`) — this cannot be done unattended.

---

## Part 5 — Order, cost, verification

| Step | Runs | CPU-h | Wall | Gate |
|---|---|---|---|---|
| 0a. Land `RecordLimitation`; run goldens to prove default-off is byte-identical | 0 | 0 | ~1 h human | goldens pass |
| 0b. Smoke the chain: `probeCover` at `NCells=30`, `tsim=0.2`, both arms, then `collectCover` | 2 | 0.05 | 3 min | collector emits pairs |
| 1. **Stage A** light ladder | 5 | 2.1 | 25 min | G-A1/A2/A3 |
| 2. **X-BIO** (+ DOM calibration runs) | 3+3 | 0.6 | 25 min | ≥ 1.5 log spread |
| 3. **X-TSPLIT** (+ winter mature cache) | 4+1 | 0.9 | 35 min | Δ_T(on) ≥ 5×Δ_T(off) |
| 4. **Stage B** nutrient × influent | 14 | 5.8 | 25–50 min | G-B1…B4 |
| 5. **Stage B′** if and only if G-B4 | 4 | ~1.5 | 30 min | — |
| 6. **Stage C** at 104 d | 4 | 5.5 | 85 min | joint criterion |
| 7. Collect, tabulate, figures | — | ~0 | 30 min | — |

**Total ≈ 17–19 CPU-hours, ≈ 4 h cluster wall across dependent arrays, ≈ 3 h human gate review.**

### Order rationale

1. **Diagnostic before sweep.** Landing `RecordLimitation` costs no compute and changes what
   Stage B should be. If the stack turns out N-limited rather than P-limited, a P-only sweep
   relieves the wrong nutrient and returns a confident null for the wrong reason.
2. **Stage A before anything structural.** There is a real chance the 4 % is a metric artefact
   of integrating φ_b instead of top-layer biomass carbon. Twenty-five minutes to rule it in or
   out, and it would be bad to change model structure without checking.
3. **X-BIO and X-TSPLIT can run concurrently with Stage A** — different caches, no dependency.
   They are also the two that directly answer a referee, so they should not queue behind the
   diagnostic work.
4. **Stage B′ is last resort, and its null is publishable.** If A and B are both flat, the
   honest finding is that the model's euphotic depth (~5 mm) is below one grid cell, so it
   structurally cannot grow a schmutzdecke and the covered/uncovered contrast is mediated
   entirely by supernatant settling. That is a real structural conclusion for the limitations
   section — B′'s two counterfactuals are what make it defensible rather than speculative.

### Verification

- **Regression:** `Pkg.test`-equivalent MATLAB goldens pass unchanged with `RecordLimitation`
  off; `rx` bit-identical before/after the `evaluateReactions` change.
- **Self-control:** Stage A's S = 0.01 cell at the manuscript stack must recover the known
  ~4 % / 0.97 ratio. If it does not, the harness is wrong, not the model.
- **Registry integrity:** `collectCover`'s `isequal(cfg, cases(i))` assertion fires on any
  stale artefact.
- **Mass guard:** per-arm ε-weighted residuals within 5 % of each other, else the pair is
  rejected.
- **End-to-end:** `collectCover("A")` produces `pairs.csv` with one row per `pairKey` and a
  populated `eps_L` column.

---

## Amendment 2026-08-22 — Part 3 executed; Stage B respecified, Stage B′ promoted

`RecordLimitation` landed and was run on both stacks
(`.claude/decisions/2026-08-22-liebig-limitation-diagnostic.md`). Three changes follow.

1. **The nutrient prior was wrong twice over.** Manuscript stack is NH₄-limited; corrected
   stack is **IC**-limited (supernatant biofilm) and **HPO₄**-limited (flowing). Neither is
   the phosphorus-recycling story assumed in `2026-08-20-winter-summer-flip.md`. Stage B's
   `+N` arm relieves nothing on the corrected stack and is replaced by **`+IC`**; add a
   **`+ALL`** saturating arm to establish the ceiling.

2. **Euphotic depth is measured, not inferred.** Î falls from 0.0608 at z = 0 to
   **1.39e-27** one cell into the bed, on both stacks. Bed photosynthesis is absent, not
   small. Stage B′ moves from "conditional on G-B4" to **scheduled**, and runs in parallel
   with Stage B rather than after it.

3. **The two blockers are regionally separated, so B and B′ are complementary.** Corrected
   stack, top 0–2 cm: min-Monod **0.703** — nutrients are essentially unlimited there, and
   light is 1e-27. Supernatant: nutrients throttle to ≈ 0.005, and light is present. So
   covering fails in the bed for lack of light and in the supernatant for lack of nutrients.
   Report both limbs; do not expect one factor to carry the result.

Revised gate logic: G-B4 ("nothing moves → euphotic depth gates") is superseded — B′ runs
regardless. G-B1/B2/B3 still apply to the supernatant limb.
