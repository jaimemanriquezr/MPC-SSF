# EXPERIMENTS.md — reproducible record of accepted runs

One section per experiment: everything between the launcher and the result. Rules, so the
record stays reproducible without reading the transcript:

* **Every run goes here** — accepted, superseded, aborted or null. A null result is a result;
  a run that clogged tells the next person not to repeat it. Interpretation and the reasons
  behind a choice stay in `.claude/JOURNAL.md` and `.claude/decisions/`; the *facts* live here.
* Each entry records, at minimum: **date/time logged**, the **verbatim command(s)**, every
  parameter that differs from the preset *and* the preset values in force at the time (they
  change — cite the decision file), **N**, run length and leg structure, **wall time per arm**,
  the **output files**, the **figures produced** (path, and what each shows), the measured
  outcome, and the caveats that must travel with the numbers.
* Overrides passed as `NaN` in a `rec` field mean "preset value"; write the number the preset
  resolved to, not `NaN`.
* If a run produced no figure, say so explicitly ("figures: none").

---

## E1 — Field and 2× field influent at the working set (2026-08-26)

**Logged** 2026-08-26 16:24 (local). **Status:** best set to date on the literature targets;
supersedes the ζ₀ = 5 / √-detachment sets of the same day.

### Commands

```matlab
% 2x field, lit / dark  (E1a / E1b)
probeChain('fld2x_lit',  Zeta0=1, Zeta1=0.27, NCells=500, Days=10, MaxDt=5e-5, LegDays=10, ...
           FramesPerLeg=241, LightScale=1, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, ...
           DetachForm="linear", Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])
probeChain('fld2x_dark', ... LightScale=0 ... )   % otherwise identical
% field, lit / dark  (E1c / E1d): Influent=[1.5e-4, 5.0e-4, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3]
probeChain('fld1x_lit', ...);  probeChain('fld1x_dark', ...)
```

### Parameters

| group | parameter | value | provenance |
|---|---|---|---|
| **influent** (kg/m³) | HET / PHO | **3.0e-4 / 1.0e-3** (2× field) · **1.5e-4 / 5.0e-4** (field) | field = Campos2006b Fig. 1(a), Kempton Park influent Chl-a 3–7 µg/L at C:Chl-a 50 |
| | POM / PAT | 0 / 0 | manuscript |
| | O₂ / IC / NH₄ / HPO₄ / DOM | 9.10e-3 / 6.23e-3 / 2.0e-5 / **5.0e-6** / **1.0e-3** | saturation; manuscript; Chan2018 (<10 µg/L total P); DOM = biodegradable fraction of a 3–4 mg/L DOC |
| **cohesion** | ζ₀ / ζ₁ / κ | **1** / 0.27 / 1e-6 | ζ₀ = 1 chosen to keep biomass out of the supernatant (Jaime, 2026-08-26); ζ₁ from the tumour-potential mapping |
| **detachment** | form / coefficient / scale | **linear**, 0.14·\|v_f\|/18 per day, ×1 | equal to the manuscript √-form at 18 m/d; linear needed to cap the z = 0 pile-up at ζ₀ = 1 |
| **transfer** | flowing↔enclosed liquids | constant form, **×10** (600 → 6000 /d; DOM 300 → 3000) | uncited multiplier; equivalent to b = D/L_f² with L_f ≈ 170 µm (Wolf2007 D_O₂ = 1.73e-4 m²/d) |
| **kinetics** (preset `modelLund`, 2026-08-25 audit) | μ_HET / θ | 2.00 /d · 1.0725 | Reichert2001 k_gro,H,aer; θ = e^β_H |
| | μ_PHO / θ | 2.00 /d · 1.047 | Reichert k_gro,ALG; θ = e^β_ALG |
| | d_HET / θ | 0.40 /d · 1.0725 | Wolf2007 b_ina,H |
| | d_PHO / θ | 0.276 /d · 1.080 | Campos2006 k_ra (respiration + excretion); single loss term, **no respiration reaction** |
| | k_hyd / θ / K_hyd | 3.00 /d · 1.0725 · 0.1 | Wolf/Reichert; Wolf K_S,h,X |
| | K_HET: O₂ / DOM / NH₄ / HPO₄ | 2.0e-4 / **3.0e-4 (override)** / 1.0e-6 / 2.0e-5 | Reichert; K_DOM override (preset 4.0e-3 = Wolf ASM value, too high for SSF) |
| | K_PHO: IC / NH₄ / HPO₄ | 1.2e-3 / 2.0e-5 / **1.0e-6 (override)** | Wolf; K_HPO₄ override = Campos2006 ksp low end (preset 2.0e-5 = Reichert, river model) |
| | dark-growth floor | 0 | artefact term retired 2026-08-26 |
| **light** | irradiance | 0.8·max(sin 2π(t − 13/48) + 0.62, 0)/1.62, normalised (I_opt = 1), Steele | manuscript curve; `LightScale` 1 (uncovered) / 0 (covered) |
| | η_water / η_sand | 1 /m · **1500 /m** | Gallegos2000 / Kühl1994 |
| **filter** | depth / supernatant / q / T | 1 m / 1 m / 7.2 m/d / 19 °C | manuscript (Schijven2013 pilot) |
| | ε₀ / δ (SandRoughness) | 0.4 / **5.0e-3 m** | manuscript |
| | attachment (sand / matrix) | 547 /d each (POM 0) | manuscript |
| | β (BiofilmPorosity) / osmosis | 0.99 / 1.0e-7 | Melo2005 / manuscript |
| | ρ_P / ρ_L | 1117 / 998 kg/m³ | Lewis2014 / water — **wet-mass convention, unresolved (see caveats)** |
| **numerics** | N / dz | **500** / 2.0 mm | N = 200 has no cell in the roughness layer and is blind to light |
| | T / legs / frames | 10 d / 1 × 10 d / 241 (hourly) | |
| | scheme | shin, Neumann cohesion BC, implicit dispersion, implicit osmosis, upwinded | `SolverOptions` defaults |
| | AdaptiveMaxDt / CFL | 5e-5 d / 0.99 | CFL-bound at N ≥ 500 |

### Wall time (MacBook, single core per arm, 4 arms concurrent + 3 other jobs)

| run | tag | wall (s) | s per simulated day |
|---|---|---|---|
| E1a 2× field lit | `fld2x_lit` | 941 | 94 |
| E1b 2× field dark | `fld2x_dark` | 860 | 86 |
| E1c field lit | `fld1x_lit` | 946 | 95 |
| E1d field dark | `fld1x_dark` | 866 | 87 |

### Output and figures

Data: `analysis/probes/data/chain/chain_{fld1x,fld2x}_{lit,dark}_leg1.mat`
(v7.3: `results` — the full `Results` object incl. `SolverOptions`; `results_py` — flattened
arrays for Python; `rec` — metrics and profile history; `snap` — restart state).

| figure | shows |
|---|---|
| `analysis/results/figures/field_load_10d.png` | six panels: (a) effluent O₂ with Elemo's window, both loads lit/dark plus the 3e-3 set for scale; (b) BDOC removal with the Campos band; (c) H/H₀; (d) φ_b near the surface at 10 d; (e) phototrophs in the roughness layer, lit vs dark; (f) O₂, DOC and H/H₀ against load at 10 d |

### Outcome at 10 d

| observable | field (5e-4) | 2× field (1e-3) | literature target |
|---|---|---|---|
| O₂ consumed (mg/L) | 1.18 | **1.97** | 2–5, never < 3 (Elemo2024, mature filter) |
| effluent O₂ / minimum | 7.92 / 7.85 | 7.13 / 7.03 | filtrate 5–8 |
| BDOC removal | 58 % | **78 %** (≈ 23 % of total DOC) | 23–25 % of total DOC (Campos2002) |
| particulate removal | 98.9 % | **99.0 %** | 98.7–99.3 % (Bae2023) |
| H/H₀ (bed) | 1.03 | 1.06 | 1.3–1.8 at 20 d → 2–4× at 55 d (Demir2017) — **not met** |
| top-2 cm biomass | 72 µg C/g dry (18 wet) | 140 (35 wet) | 60 (Campos2002, d97) |
| biomass maximum | at the sand surface | at the sand surface | top 1.5–2 cm |
| supernatant biomass | 0.0 % | 0.0 % | thin schmutzdecke only |
| lit/dark: bed · top-2 cm · PHO above sand | 1.006 · 1.004 · 1.31× | 1.012 · 1.010 · **1.42×** | 4× in the top 2 cm (Campos2002) — **not met** |

### Caveats (must accompany any use of these numbers)

1. **Wet vs dry mass is unresolved.** The state variable is wet biomass (ρ_P = 1117 from a
   cell-volume × density derivation) while the RWQM1 stoichiometry is per dry organic matter;
   yields are ~4× too strong per kg of state variable and the µg C/g figures move by the same
   factor. Under the wet reading, 1e-3 *is* Campos's field influent in our units.
2. **10 d is a transient.** Effluent O₂ is still falling in all four arms; Elemo's filter is
   mature (15 months, post-skimming steady state reached after 14 d). These values are lower
   bounds on the eventual demand.
3. **Head loss is far below Demir's**, and no arm has been run to her 55 d horizon.
4. **The light contrast is real but small** and lives in the 4 mm roughness layer; the model
   has no mechanism retaining phototrophs in the light (see `.claude/CRITIQUE.md` §6).
5. K_DOM 3e-4 and the ×10 transfer are the two uncited numbers in the set.

---

## E2 — Heterotroph growth-rate ladder (2026-08-26)

**Logged** 2026-08-26 16:05. **Status:** rejected as a lever — μ_HET trades O₂ against DOC
removal one-for-one. Recorded because it closes the question.

### Commands

```matlab
% working set, 5 d, lit; only MuHET differs (NaN = preset 2.0)
probeChain('muhet_2p0',  Zeta0=1, Zeta1=0.27, NCells=500, Days=5, MaxDt=5e-5, LegDays=5, ...
           FramesPerLeg=121, LightScale=1, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, ...
           DetachForm="linear", MuHET=NaN, ...
           Influent=[8.04e-4, 3.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])
probeChain('muhet_1p0',  ... MuHET=1.0  ... );  probeChain('muhet_0p5', ... MuHET=0.5 ... )
probeChain('muhet_0p25', ... MuHET=0.25 ... )
```

### Parameters

As E1 (same preset, ζ₀ = 1, linear detachment ×1, transfer ×10, K_DOM 3e-4, K_HPO₄ 1e-6,
N = 500, δ 5 mm, η_sand 1500, lit) with two differences: **influent PHO/HET = 3.0e-3 / 8.04e-4**
(the ladder rung between field and manuscript, not a cited value) and **μ_HET ∈ {2.0 (preset,
Reichert), 1.0 (Jaime's 08-20 audit, 0.042 × 24), 0.5, 0.25}**; Campos2006's bacterial range
is 0.70–1.0, so 0.5 and 0.25 are below every cited value. Run length **5 d**, one leg,
121 frames (2-hourly).

### Wall time

| tag | μ_HET | wall (s) |
|---|---|---|
| `muhet_2p0` | 2.0 | ≈ 640 |
| `muhet_1p0` | 1.0 | ≈ 640 |
| `muhet_0p5` | 0.5 | ≈ 640 |
| `muhet_0p25` | 0.25 | ≈ 640 |

(all four launched together and completed within the same ~11 min window; per-arm values are
in each file's `rec.wall`.)

### Output and figures

Data: `analysis/probes/data/chain/chain_muhet_{2p0,1p0,0p5,0p25}_leg1.mat`.

| figure | shows |
|---|---|
| `analysis/results/figures/muhet_ladder_5d.png` | (a) effluent O₂ per rung; (b) BDOC removal — negative below μ = 2.0; (c) bed biofilm; (d) O₂ down the bed at 5 d; (e) heterotroph profiles; (f) the O₂-vs-DOC trade-off against μ_HET |

### Outcome at 5 d

| μ_HET | O₂ consumed | BDOC removal | bed | HET share |
|---|---|---|---|---|
| 2.0 | 2.38 mg/L | +83 % | 0.00922 | 40 % |
| 1.0 | 1.22 | −20 % | 0.00812 | 32 % |
| 0.5 | 0.17 | −123 % | 0.00729 | 24 % |
| 0.25 | −0.21 (net producer) | −158 % | 0.00696 | 20 % |

**Conclusion:** lowering μ_HET improves the oxygen balance only by leaving captured carbon
unmineralised — the filter becomes a DOM source. Only μ_HET = 2.0 puts both O₂ and DOC in
their literature bands. The O₂ demand is set by the organic load, not the heterotroph rate;
the load is the lever (see E1).

---

## E3 — ζ₀ = 1 with linear detachment (2026-08-26)

**Logged** 2026-08-26 15:30. **Status:** accepted — established the working set that E1/E2 use.
Supersedes ζ₀ = 5 with √-detachment.

### Commands

```matlab
probeChain('z01_set_lin_lit',  Zeta0=1, Zeta1=0.27, NCells=500, Days=20, MaxDt=5e-5, LegDays=10, ...
           FramesPerLeg=241, LightScale=1, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, ...
           DetachForm="linear", ...
           Influent=[8.04e-4, 3.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])
probeChain('z01_set_lin_dark', ... LightScale=0 ... )
```

### Parameters

As E1 with PHO/HET = 3.0e-3 / 8.04e-4 and **20 d in two 10 d legs**. Previous ζ₀ = 1 runs
with the √-detachment form clogged at 14–17 d (`chain_z0_1_n500_hs_{neu,dir}_leg2_ABORTED.mat`);
the linear form caps the z = 0 pile-up (k_det 0.9–2.8 /d there against 0.35–0.65 with √).

### Wall time

| tag | legs | wall (s) |
|---|---|---|
| `z01_set_lin_lit` | 892 + 1115 | 2007 total (100 s per simulated day) |
| `z01_set_lin_dark` | 817 + 1105 | 1922 |

### Output and figures

Data: `analysis/probes/data/chain/chain_z01_set_lin_{lit,dark}_leg{1,2}.mat`.

| figure | shows |
|---|---|
| `analysis/results/figures/zeta0_1_linear_lit_vs_dark.png` | (a) bed and supernatant biofilm; (b) the surface cell — no clog; (c) effluent O₂, days 16–20; (d) φ_b near the surface; (e) φ_b full bed; (f) phototroph and heterotroph profiles |

### Outcome at 20 d

| | lit | dark |
|---|---|---|
| bed / supernatant share | 0.0181 m / **0.3 %** | 0.0173 / 0.0 % |
| φ_b(0) (bounded, no clog) | 0.514 | 0.453 |
| H/H₀ | 1.29 | 1.28 |
| O₂ mean / minimum | 3.88 / 3.58 | 3.80 / 3.75 |
| top-2 cm | 352 µg C/g dry | 334 |
| lit/dark: bed · top-2 cm · PHO above sand | 1.045 · 1.053 · 2.46× | — |

**Note:** the interface cell is still a local maximum (0.51 against neighbours at 0.38/0.32) —
a bounded version of the ζ₀ = 1 pile-up. Whether the linear form holds it at 40–60 d is untested.
The linear detachment form has no decision file yet; the physical argument is shear ∝ v/r
through a closing pore.

---

## E4 — 2× field influent to 60 days (2026-08-26)

**Logged** 2026-08-26 19:05. **Status:** accepted; the first run of this project to reach a
steady state and to cover Demir2017's 55-day horizon. Extends E1a/E1b (same tags, legs 2–6
appended to the existing 10 d leg 1).

### Commands

```matlab
probeChain('fld2x_lit',  Zeta0=1, Zeta1=0.27, NCells=500, Days=60, MaxDt=5e-5, LegDays=10, ...
           FramesPerLeg=241, LightScale=1, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, ...
           DetachForm="linear", ...
           Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])
probeChain('fld2x_dark', ... LightScale=0 ... )   % otherwise identical
```

Leg 1 (0–10 d) was already on disk from E1 and was detected and skipped; legs 2–6 ran.

### Parameters

Identical to E1a/E1b (see that entry for the full table): 2× field influent
(HET 3.0e-4 / PHO 1.0e-3 kg/m³), DOM_in 1.0e-3, HPO₄_in 5.0e-6; ζ₀ = 1, ζ₁ = 0.27, κ = 1e-6;
**linear detachment** 0.14·|v_f|/18, ×1; constant liquid transfer ×10; K_DOM 3.0e-4,
K_HPO₄,PHO 1.0e-6; preset kinetics of the 2026-08-25 audit (μ_HET 2.0/θ 1.0725, μ_PHO 2.0/θ 1.047,
d_HET 0.40, **d_PHO 0.276/θ 1.08, no respiration reaction**, k_hyd 3.0); dark floor 0;
η_sand 1500 /m, δ 5 mm; N = 500, MaxDt 5e-5, shin scheme, Neumann cohesion BC.
Run length **60 d in six 10 d legs**, 241 hourly frames per leg.

### Wall time

| leg | t (d) | lit (s) | dark (s) |
|---|---|---|---|
| 1 | 0–10 | 940 | 859 |
| 2 | 10–20 | 907 | 899 |
| 3 | 20–30 | 1077 | 1098 |
| 4 | 30–40 | 1095 | 1121 |
| 5 | 40–50 | 1075 | 1110 |
| 6 | 50–60 | 1054 | 1078 |
| **total** | | **1.71 h** | **1.71 h** |

≈ 105 s per simulated day at N = 500; cost per day rises ~15 % from leg 2 to leg 3 and then
flattens with the biomass.

### Output and figures

Data: `analysis/probes/data/chain/chain_fld2x_{lit,dark}_leg{1..6}.mat` — a checkpoint every
10 d, each with `results`, `results_py`, `rec` (241 hourly frames) and `snap`.

| figure | shows |
|---|---|
| `analysis/results/figures/fld2x_60d.png` | (a) bed biofilm, steady from ~30 d; (b) effluent O₂ with Elemo's window; (c) H/H₀ with Demir's 55 d band and the 55 d marker; (d) BDOC removal; (e) φ_b profiles at 10–60 d; (f) top-2 cm biomass and the lit/dark phototroph ratio |
| `analysis/results/figures/fld2x_60d_partial.png` | the same at 40 d (superseded, kept for the record) |

### Outcome

| t (d) | bed lit / dark (ratio) | H/H₀ | O₂ mean / min | consumed | BDOC | top-2 cm (wet) | PHO above sand |
|---|---|---|---|---|---|---|---|
| 10 | 0.00548 / 0.00542 (1.012) | 1.06 | 7.13 / 7.03 | 1.97 | 78 % | 35 | 1.42× |
| 20 | 0.00732 / 0.00714 (1.025) | 1.09 | 6.39 / 6.32 | 2.71 | 87 % | 43 | 2.05× |
| 30 | 0.00780 / 0.00754 (1.034) | 1.09 | 6.20 / 6.13 | 2.90 | 88 % | 44 | 2.36× |
| 40 | 0.00792 / 0.00764 (1.038) | 1.09 | 6.15 / 6.08 | 2.95 | 89 % | 44 | 2.48× |
| 50 | 0.00796 / 0.00766 (1.040) | 1.10 | 6.14 / 6.07 | 2.96 | 89 % | 45 | 2.52× |
| 60 | 0.00797 / 0.00766 (1.040) | 1.10 | 6.14 / 6.07 | 2.96 | 89 % | 45 | 2.54× |

**Steady state reached:** over the last 10 days bed biomass changes +0.14 % and effluent O₂
by −0.006 g/m³. Supernatant biomass stays at 0.0 % of the total throughout.

**Against the targets:** O₂ consumption **2.96 mg/L**, inside Elemo2024's 2–5 window with the
minimum at 6.07 (floor 3.0) — *met, at steady state, against a mature-filter measurement*.
BDOC removal 89 % ≈ 27 % of a 3.3 mg/L total DOC — *met* (Campos2002 23–25 %), marginally high.
Bed biomass and head loss monotone throughout — *met*. Biomass maximum at the sand surface,
supernatant empty — *met*. **H/H₀ plateaus at 1.10 against Demir's 2.1–4.2 at 55 d — not met,
and now measured rather than extrapolated.** Top-2 cm 45 µg C/g wet vs Campos's 60 at d97 —
close on the wet reading, 4× high on the dry one.

**Light contrast grows while everything else saturates:** phototrophs above the sand reach
**2.54×** dark and were still rising at 60 d, and the bed ratio widens 1.012 → 1.040. This is
the signature of a *grown* surface layer rather than trapping, which would have saturated with
the rest — the clearest evidence so far that the light response is photosynthetic.

### Caveats

1. The set does not clog. At this load the filter reaches a steady state at H/H₀ = 1.10 and
   stays there; no filter run ends. Reproducing Demir's or Campos's clogging needs a higher
   particulate load (3e-3 gave 1.29 at 20 d, still short) or the missing surface-mat
   mechanisms (`.claude/CRITIQUE.md` §6, §12).
2. Wet/dry mass convention still unresolved; the µg C/g figures move by ~4×.
3. K_DOM 3e-4 and the ×10 transfer remain the two uncited numbers.
4. The lit/dark contrast is 4 % in the bed — real and growing, but far from Campos's 4×.

---

## E5 — Phototroph growth rate ×10 and ×20 at N = 500 (2026-08-26)

**Logged** 2026-08-26 19:20. **Status:** diagnostic, negative result. μ_PHO is not the lock on
the light contrast. Both rates are far outside any cited value (Reichert k_gro,ALG = 2.0,
Campos 1.0–3.0) and are **not** candidate settings.

### Commands

```matlab
probeChain('mu500_x10_lit',  Zeta0=5, Zeta1=0.27, NCells=500, Days=20, MaxDt=5e-5, LegDays=10, ...
           FramesPerLeg=241, LightScale=1, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, MuPHO=20.0, ...
           Influent=[2.68e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])
probeChain('mu500_x10_dark', ... LightScale=0 ...);  probeChain('mu500_x20_lit', ... MuPHO=40.0 ...)
probeChain('mu500_x20_dark', ... MuPHO=40.0, LightScale=0 ... )
```

### Parameters

As E1 **except**: ζ₀ = **5** (these were launched before the ζ₀ = 1 decision of E3), √-detachment
(the preset form, also pre-E3), influent HET 2.68e-4 / PHO 1.0e-3, and **μ_PHO = 20 or 40 /d**
(preset 2.0). N = 500, 20 d in two 10 d legs, 241 hourly frames. K_DOM 3e-4, K_HPO₄ 1e-6,
transfer ×10, η_sand 1500, δ 5 mm, d_PHO 0.276, no respiration term.

### Wall time — the cost of a stiff reaction term

| tag | μ_PHO | leg 1 (s) | leg 2 (s) | total |
|---|---|---|---|---|
| `mu500_x10_lit` | 20 | 5035 | 5066 | **2.81 h** |
| `mu500_x10_dark` | 20 | 2019 | 1576 | 1.00 h |
| `mu500_x20_lit` | 40 | 9371 | 8708 | **5.02 h** |
| `mu500_x20_dark` | 40 | 4118 | 2475 | 1.83 h |

Against ~0.55 h per 20 d at the preset μ_PHO = 2: a 10× rate costs 5× the wall time and a 20×
rate 9×, because the reaction term drives the adaptive step down. The lit arms cost 2.8× their
dark twins — the stiffness is in the growth term, which only acts where there is light.

### Output and figures

Data: `analysis/probes/data/chain/chain_mu500_x{10,20}_{lit,dark}_leg{1,2}.mat`.
**Figures: none** (the result is a table; no figure was produced).

### Outcome at 20 d

| μ_PHO | bed lit / dark | lit/dark bed | supernatant lit/dark | PHO above sand | H/H₀ | O₂ |
|---|---|---|---|---|---|---|
| 2 (preset, E4 at the same influent) | 0.00732 / 0.00714 | 1.025 | — (0.0 %) | 2.05× | 1.09 | 6.39 |
| 20 (×10) | 0.00749 / 0.00697 | 1.075 | 1.68× | 3.12× | 1.09 | 6.49 |
| 40 (×20) | 0.00755 / 0.00697 | 1.083 | 1.64× | 3.14× | 1.09 | 6.50 |

**A 20× increase in the phototroph growth rate buys a lit/dark bed ratio of 1.08 instead of
1.03, and saturates**: 40 /d is indistinguishable from 20 /d in every column. Supernatant
biomass stays at 0.5 % of the total. Head loss and effluent O₂ do not move at all.

**Conclusion.** The light contrast is limited by *where phototrophs can live*, not by how fast
they grow: with only ~4 mm of lit, attaching roughness layer, even an absurd growth rate cannot
build a mat or move the bed. This is the same conclusion as the K_HPO₄ ladder and the δ = 2 cm
arm, from the opposite direction, and it closes μ_PHO as a lever. See `.claude/CRITIQUE.md`
§6 (no retention of phototrophs in the lit zone) and §12 (no buoyant mat export).

---

## E6 — μ_PHO ×10 on the working set (ζ₀ = 1, linear detachment) — COSMOS (2026-08-27)

**Logged** 2026-08-27. **Status:** accepted; supersedes E5's ζ₀ = 5 arms for comparison with E4.
**Ran on LUNARC cosmos**, job array `3544074` (partition lu48, node cn063, account lu2026-2-100).

### Commands

```bash
# cosmos: MPC-SSF-manuscript/slurm/mu_z01.sbatch  (array 0-1; 0 = lit, 1 = dark)
sbatch slurm/mu_z01.sbatch
```
```matlab
probeChain("mu_z01_x10_lit", Zeta0=1, Zeta1=0.27, NCells=500, Days=20, MaxDt=5e-5, LegDays=10, ...
           FramesPerLeg=241, LightScale=1, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, ...
           DetachForm="linear", MuPHO=20.0, ...
           Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])
```

### Parameters

Identical to E4 (2× field influent, ζ₀ = 1, linear detachment ×1, transfer ×10, K_DOM 3e-4,
K_HPO₄ 1e-6, N = 500, δ 5 mm, η_sand 1500, audited preset, no respiration reaction) with the
single change **μ_PHO = 20 /d (×10 the preset 2.0)**. 20 d in two 10 d legs, 241 hourly frames.
μ_PHO = 40 was not repeated: E5 showed it identical to 20 in every column.

**Prerequisite:** the cosmos tree was two days stale (no `TransferForm`/`DetachForm`/`MuHET`
options; pre-audit `modelLund`). `src/` and `analysis/` were rsynced and the corrected preset
verified remotely (`d_PHO = 0.276`, `k_hyd = 3.00`, `K_NH₄,PHO = 2.0e-5`) before submitting.

### Wall time

| task | arm | leg 1 (s) | leg 2 (s) | total |
|---|---|---|---|---|
| `3544074_0` | lit | 8889 | 7485 | **4:33:25** |
| `3544074_1` | dark | — (resumed) | 3694 | **1:02:02** |

Cosmos is ~1.8× slower per leg than the local machine for this arm (8889 s vs 4371 s locally),
so cosmos is worth it for parallel breadth and for jobs that outlive the session, not for speed.

### Output and figures

Data: `analysis/probes/data/chain/chain_mu_z01_x10_{lit,dark}_leg{1,2}.mat`, fetched by rsync
from `cosmos:MPC-SSF-manuscript/analysis/probes/data/chain/`. **Figures: none.**

### Outcome

| t (d) | bed lit / dark | ratio | top-2 cm lit/dark | PHO above sand | H/H₀ | O₂ lit |
|---|---|---|---|---|---|---|
| 10 | 0.00561 / 0.00542 | 1.036 | 1.026 | 1.92× | 1.07 | 7.21 |
| 20 | 0.00765 / 0.00714 | **1.071** | **1.093** | **2.74×** | 1.09 | 6.43 |
| *E4 reference, μ_PHO = 2* | *0.00732 / 0.00714* | *1.025* | *1.010* | *2.05×* | *1.09* | *6.39* |

A tenfold phototroph growth rate lifts the bed contrast from 1.025 to **1.071** and the top-2 cm
from 1.010 to **1.093** — roughly a threefold gain in the contrast for a tenfold rate, and it
still leaves the supernatant empty (0.000000 m at 20 d, against 0.5 % at ζ₀ = 5 in E5) and
head loss and effluent O₂ unchanged. Consistent with E5's conclusion: growth rate modulates the
contrast weakly and saturates; the binding constraint is the ~4 mm of lit, attaching zone.

### Caveat — data lost by operator error

Local copies of `chain_mu_z01_x10_*` were deleted with `rm -f` on the assumption they were
partial, when the dark arm had completed both legs and the lit arm leg 1 (~2.5 h of compute).
The cosmos run regenerated them. Verify before deleting.

---

## E8 — Load ladder at the working set: 1× and 4× field influent, 60 d — COSMOS (2026-08-27)

**Logged** 2026-08-27 21:45 (launch) / 23:0x (outcome, below). **Status:** running at the time
of logging; partial results only. **Ran on LUNARC cosmos**, job array **`3549851`** (partition
lu48, account lu2026-2-100, nodes cn096/cn159), submitted 21:24 CEST.

Purpose: E4/E7 established that the 2× field load reaches a steady state at H/H₀ = 1.10 and
never clogs. The ladder asks whether the load is the lever: 1× field and 4× field, lit and dark,
otherwise byte-identical to E4. The 2× rung is E4/E7 itself and is not repeated.

### Commands

```bash
# cosmos: MPC-SSF-manuscript/slurm/ladder_n500.sbatch  (array 0-3)
sbatch slurm/ladder_n500.sbatch
```
```matlab
% task 0/1: 1x field (HET 1.5e-4, PHO 5.0e-4); task 2/3: 4x field (6.0e-4, 2.0e-3)
probeChain("fld1x_lit", Zeta0=1, Zeta1=0.27, NCells=500, Days=60, MaxDt=5e-5, LegDays=10, ...
           FramesPerLeg=241, LightScale=1, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, ...
           DetachForm="linear", ...
           Influent=[1.5e-4, 5.0e-4, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])
% fld1x_dark  LightScale=0;  fld4x_lit / fld4x_dark  Influent(1:2) = [6.0e-4, 2.0e-3]
```

### Parameters

Identical to E4 in every option except `Influent(1:2)` and `LightScale` — see the E4 entry for
the full table (ζ₀ = 1, ζ₁ = 0.27, κ = 1e-6, linear detachment ×1, transfer ×10, K_DOM 3.0e-4,
K_HPO₄ 1.0e-6, audited preset kinetics with no respiration reaction, dark floor 0, η_sand 1500,
δ 5 mm, N = 500, MaxDt 5e-5, shin scheme, Neumann cohesion BC, 19 °C, 60 d in six 10 d legs,
241 hourly frames per leg).

**The `fld1x` arms are a RESUME**, not a fresh chain: `chain_fld1x_{lit,dark}_leg1.mat` (0–10 d)
already existed from E1c/E1d with exactly these options, so those arms skipped leg 1 and
computed legs 2–6. Confirmed in `slurm/logs/ladder_3549851_0.out`
("leg 1 already done (t = 10 d) -- skipped"). The `fld4x` arms started from t = 0.

### Wall time

Not available for any arm at logging time (the array was still running). Estimated ~190 s per
simulated day on cosmos (E6 measured 1.8× the local 105 s/d), i.e. ~32 min per 10 d leg; three
of the four arms share node cn096 with each other and with Agent A's pulse array, so per-leg
times will be worse than that.

### Output and figures

Data: `analysis/probes/data/chain/chain_{fld1x,fld4x}_{lit,dark}_leg*.mat` on
`cosmos:MPC-SSF-manuscript/`, fetched by rsync as legs land. **Figures: none.**
Launcher: `slurm/ladder_n500.sbatch`; logs `slurm/logs/ladder_3549851_{0..3}.{out,err}`.

### Outcome

Scored with `analysis/scoreRun.m` at whatever time each arm had reached. **See the table
appended at the end of this entry**; any arm listed as "not finished" is exactly that — the
number given is at the reached time, not at 60 d.

**The arms had NOT finished when this entry was written.** Reached times below are what was on
disk at 22:15 CEST; the 60 d target was not met by any arm inside the session. Per-leg wall time
on cosmos: fld4x leg 1 1273 s (lit) / 1204 s (dark), fld1x leg 2 1620 s / 1609 s — i.e. ~21–27
min per 10 d leg, so the remaining legs land at roughly half-hour intervals and the job (12 h
limit) will complete on its own.

| arm | reached t (d) | bed ∫εφ_b (m) | total column biomass (kg/m²) | H/H₀ | O₂ consumed (mg/L) | O₂ min | top-2 cm wet / dry (µg C/g) | PHO above sand lit/dark | **score** |
|---|---|---|---|---|---|---|---|---|---|
| `fld1x_lit` | 20 | 0.00426 | 0.04755 | 1.047 | 1.83 | 7.23 | 23 / 93 | 1.76× | **7/10** (3/2/2) |
| `fld1x_dark` | 20 | 0.00419 | 0.04679 | 1.047 | 1.84 | 7.24 | 23 / 91 | — | **6/10** (3/1/2) |
| *`fld2x_lit` (E4/E9, for scale)* | *104* | *0.00798* | *0.08909* | *1.095* | *2.97* | *6.06* | *42 / 167* | *2.55×* | ***10/10*** |
| `fld4x_lit` | 20 | 0.01300 | 0.14527 | 1.181 | 4.11 | 4.82 | 65 / 260 | 2.36× | **9/10** (3/3/3) |
| `fld4x_dark` | 20 | 0.01250 | 0.13959 | 1.175 | 4.16 | 4.92 | 62 / 248 | — | **8/10** (3/2/3) |
| *`fld4x_lit` at 10 d* | *10* | *0.01020* | *0.11393* | *1.139* | *3.13* | *5.79* | *59 / 234* | *1.68×* | *9/10* |

Scores from `analysis/scoreRun.m` (Campos point scored on the **wet** reading; see the caveats).
Files: `analysis/results/scores/{fld1x,fld4x}_{lit,dark}.md`.

**Justification per axis, in one line each.** Profile: every arm has its φ_b maximum at the sand
surface, decays monotonically below it and leaves the supernatant empty — 3/4 for all of them,
and the fourth point is unavailable because none has yet run 10 d past a non-zero start (the
steadiness window would reach back to the clean column). Biomass: the whole-column mass is
monotone in every arm — **an axis that cannot fail while a filter is simply ripening, so it
carries no information in this ladder** — and the Campos point separates the arms properly
(1× field 23 µg C/g wet is below the 30–120 window, 2× 42 and 4× 59 are inside, 4× essentially
on Campos's 60). Oxygen: 1× consumes 1.83 mg/L (the 1–2 margin, 2/3), 2× 2.96 and 4× 3.13 are
inside Elemo's 2–5 window with minima above the 3 mg/L floor (3/3).

**What the ladder shows so far.** Load moves everything μ_PHO could not (E5/E6). At the same
age of **20 d**, 1× / 2× / 4× field give:

| at t = 20 d | 1× | 2× (E4) | 4× |
|---|---|---|---|
| bed ∫εφ_b (m) | 0.00426 | 0.00732 | **0.01300** |
| H/H₀ | 1.047 | 1.09 | **1.181** |
| O₂ consumed (mg/L) | 1.83 | 2.71 | **4.11** |
| effluent O₂ minimum | 7.23 | 6.32 | **4.82** |
| top-2 cm (µg C/g wet) | 23 | 43 | **65** |
| PHO above sand, lit/dark | 1.76× | 2.05× | **2.36×** |

Every literature axis moves the right way with load and roughly linearly in it. The 4× arm is
already past the 2× arm's *terminal* head loss (1.10 at 104 d) while still in its first
fortnight, and its top-2 cm biomass sits on Campos's 60 µg C/g. That is the trajectory
Demir's 2.1–4.2 needs — but it is a trajectory, not an arrival: the 2× arm also looked steep at
10 d (1.06) and then plateaued, and the 4× effluent minimum has already fallen to 4.8 mg/L, so
the load that clogs may be the load that goes anoxic. Only the finished 60 d arms settle it.

### Caveats

1. **Nothing here is a converged run.** The 1× arms are at 20 d and the 4× arms at 10 d against a
   2× arm that needed 30 d to reach steady state. Every number above is a trajectory point.
2. The "whole-column biomass monotone" rubric axis is uninformative for ripening runs (see above);
   the 9/10 and 8/10 of the 4× arms are carried partly by it.
3. The Campos biomass point is scored on the **wet** reading of the state variable, with
   wet = dry / 4. That factor is quoted in E1 caveat 1 and has **no independent derivation**;
   on the dry reading the 4× arm reads 234 µg C/g and would fail the same point that the 1× arm
   fails from below. f_dry decides the direction of this whole ladder.
4. The `fld1x` arms resume an E1c/E1d leg 1 that was run on 2026-08-26; the `fld4x` arms are
   fresh. Same options, but the two arms of a rung do not have bit-identical histories.
5. `src/presets/pathogenModel.m` was edited and rsynced to cosmos by the concurrent agent at
   ~21:36, twelve minutes after this array started. probeChain reads the preset once at startup,
   so the running arms are unaffected; a resubmission would pick up the new file.

### Outcome — 60 d, all four arms finished (2026-08-28 00:20)

Job 3549851 completed: task 0 fld1x_lit 2:33:37, 1 fld1x_dark 2:30:51, 2 fld4x_lit 2:53:45,
3 fld4x_dark 2:46:39. Rescored at 60 d (`analysis/results/scores/fld{1x,4x}_{lit,dark}.md`,
wet convention).

| arm | score | P/B/O | bed ∫εφ_b | H/H₀ | O₂ consumed | O₂ min | top-2 cm wet/dry | contrast |
|---|---|---|---|---|---|---|---|---|
| fld1x_lit | **9** | 4/2/3 | 0.00483 | 1.054 | 2.09 | 6.98 | 25.1 / 100.4 | 2.10× |
| fld1x_dark | **8** | 4/1/3 | 0.00471 | 1.053 | 2.09 | 7.01 | 23.9 / 95.8 | — |
| fld4x_lit | **10** | 4/3/3 | 0.01395 | 1.193 | 4.44 | 4.49 | 66.3 / 265.3 | 2.72× |
| fld4x_dark | **9** | 4/2/3 | 0.01317 | 1.184 | 4.44 | 4.66 | 62.2 / 249.0 | — |

The 20 d partial scores rose: reaching steady state gave every arm the fourth profile point
(3→4), and the 1× oxygen consumption climbed from 1.83 (1–2 band) to 2.09, into Elemo's 2–5
window. Ladder at 60 d steady state, **1× / 2× (E4) / 4×**: bed 0.00483 / 0.00797 / 0.01395;
H/H₀ 1.054 / 1.10 / 1.193; O₂ consumed 2.09 / 2.96 / 4.44; effluent minimum 6.98 / 6.07 / 4.49;
top-2 cm 25 / 42 / 66 µg C/g wet; contrast 2.10× / 2.54× / 2.72×.

**Conclusion — load is the lever μ_PHO was not, but the model still cannot clog at realistic
O₂.** Every axis rises monotonically with load, and the trajectory worry of the 20 d caveat is
resolved: the 4× arm reaches a *terminal* H/H₀ of 1.19, not a passing value. But that is a fifth
of Demir2017's 2.1 floor, and it is bought at O₂ consumption 4.44 (top of the 2–5 window) with
the effluent minimum down to 4.5 — so a load high enough to clog would leave the O₂ window
first. Clogging and realistic oxygen are not simultaneously reachable from the load axis, the
same way E4/E9 showed they are not reachable by run length. The 4× arm is the session's
best-scoring long run (10/10) because the load lifts the top-2 cm onto Campos's 60 µg C/g (66
wet) while the other targets hold — one good operating point, not a clogging filter. f_dry still
decides the biomass sentence (dry 100–265 µg C/g).

---

## E9 — 2× field influent to 104 days (E4/E7 extension) — COSMOS (2026-08-27)

**Logged** 2026-08-27 22:00. **Status:** accepted, complete. Extends E4 (same tags, legs 7–11
appended to legs 1–6). **Ran on LUNARC cosmos**, job array **`3549014`** (lu48, node cn078),
launcher `slurm/fld2x_104d.sbatch`.

Purpose: E4 reached a steady state at 60 d but could only extrapolate to Campos2002's d97
sampling and to Demir2017's clogging horizon. 104 d is Campos's Bed 9 run length, so the
top-2 cm biomass and the head-loss trajectory become like-for-like comparisons.

### Commands

```bash
sbatch slurm/fld2x_104d.sbatch          # array 0 = lit, 1 = dark
```
```matlab
probeChain("fld2x_lit", Zeta0=1, Zeta1=0.27, NCells=500, Days=104, MaxDt=5e-5, LegDays=10, ...
           FramesPerLeg=241, LightScale=1, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, ...
           DetachForm="linear", ...
           Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])
```

### Parameters

Identical to E4; only `Days` changed, 60 → 104. Legs 1–6 were on cosmos and were detected and
skipped; legs 7–11 ran (leg 11 is the 4 d remainder).

### Wall time

| task | arm | elapsed (sacct) | exit |
|---|---|---|---|
| `3549014_0` | lit | **02:21:09** | 0:0 |
| `3549014_1` | dark | **02:21:15** | 0:0 |

44 simulated days per arm, so ≈ 193 s per simulated day on cosmos — 1.8× the local 105 s/d of
E4, exactly the E6 factor.

### Output and figures

Data: `analysis/probes/data/chain/chain_fld2x_{lit,dark}_leg{7..11}.mat` (425 MB, fetched by
rsync). Scores: `analysis/results/scores/fld2x_{lit,dark}_104d.md`.
Figures: `analysis/results/figures/repro/fig{2,4,5,6,7,10}_*.{png,pdf}` (the Manriquez2026
reproduction, see `analysis/results/figures/repro/README.md`) — no dedicated E9 figure.

### Outcome

| observable | 60 d (E4) | **104 d (E9)** | target |
|---|---|---|---|
| bed ∫εφ_b, lit / dark (m) | 0.00797 / 0.00766 | **0.00798 / 0.00766** | — |
| total column biomass, lit (kg/m²) | 0.08904 | **0.08909** = 0.08621 bed + 0.00289 roughness + 0.00000 supernatant | — |
| H/H₀ | 1.10 | **1.095** | Demir2017 2.1–4.2 at 55 d — **not met** |
| effluent O₂ mean / min (mg/L) | 6.14 / 6.07 | **6.14 / 6.06** | Elemo2024 filtrate 5–8 |
| O₂ consumed (mg/L) | 2.96 | **2.96** | 2–5 — met |
| top-2 cm (µg C/g) | 45 wet | **41.7 wet / 166.8 dry** | Campos2002 **60 at d97** — met on the wet reading |
| φ_b peak, lit / dark | — | **0.3239 / 0.2711 at z = 0** | maximum at the sand surface — met |
| PHO above sand, lit/dark | 2.54× | **2.55×** | Campos 4× in the top 2 cm — not met |
| profile change over the last 10 d | 0.14 % | **0.001 %** | — |

**Score (analysis/scoreRun.m): lit 10/10 (profile 4/4, biomass 3/3, O₂ 3/3); dark 9/10**
(the dark arm cannot earn the phototroph-contrast point, by construction).
**The 10/10 is conditional on the wet reading of the biomass state variable** (wet = dry / 4,
a factor with no independent derivation — E1 caveat 1). On the dry reading the top-2 cm is
166.8 µg C/g, 2.8× Campos, the Campos point is lost and the lit arm scores 9/10.

**What 104 d settles.**

1. The steady state of E4 is real and terminal, not a slow transient: between 60 d and 104 d the
   bed integral moves by +0.1 %, the profile by 0.001 %, and effluent O₂ by 0.01 mg/L.
2. **Demir's head loss is not reachable at this load.** H/H₀ = 1.095 at 104 d against 2.1–4.2 at
   55 d. The 60 d "not met" was not a matter of run length. E8 tests whether load fixes it.
3. Campos's d97 top-2 cm biomass is matched on the wet reading (41.7 vs 60 µg C/g, same order,
   28 % low) and missed by 2.8× on the dry one. f_dry still decides which sentence is written.
4. The light contrast has saturated: phototrophs above the sand 2.54× → 2.55×, the bed ratio
   1.040 → 1.041. The "still rising at 60 d" note in E4 does not survive to 104 d.

### Caveats

1. Wet/dry mass convention still unresolved (E1 caveat 1); the µg C/g figures move by 4×.
2. `scoreRun` computes µg C/g with 0.374 kg C per kg of state variable (0.531/1.42 for
   C₅H₇O₂N, as in `analysis/probes/probeCover.m`); the E4 session used 0.4, which is why the
   E4 row says 45 and the E9 row 41.7 for what is essentially the same state. Ratio 1.0695.
3. Effluent NH₄ and HPO₄ run at **8.4× and 6.6× the influent** at steady state (0.168 vs 0.020
   and 0.0328 vs 0.0050 g/m³): the filter is a net mineraliser of N and P. Not compared with
   any measurement yet, and not visible in the E4 outcome table because only O₂ and BDOC were
   tabulated there. Flagged for the revision.
4. K_DOM 3e-4 and the ×10 transfer remain the two uncited numbers.

---

## P1 — PAT challenge pulses on the E4 60 d mature filter (audited pathogen model) — COSMOS (2026-08-27)

**Logged** 2026-08-27. **Status:** all 8 runs complete and accepted (the two BigPulse arms
finished last and were added to this entry after the first six). **Ran on LUNARC cosmos**, job array
**3549855** (partition lu48, nodes cn159/cn039, account lu2026-2-100), 7 of the 8 tasks
co-scheduled on cn159.

First use of the audited `pathogenModel`
(`.claude/decisions/2026-08-27-pathogen-model-audit.md`): MarkerGrowth Monod set moved to the
audited heterotroph row (O₂ 2.0e-4, NH₄ 1.0e-6, **HPO₄ 2.0e-5 added**, DOM 4.0e-3) and ζ₀ /
the detachment law made options so the preset runs on the working set unchanged.

### Commands

```bash
# cosmos: MPC-SSF-manuscript/slurm/pulse_e4.sbatch   (array 0-7)
sbatch slurm/pulse_e4.sbatch
```
```matlab
probePulse("e4_lit_nom", Snapshot="chain_fld2x_lit_leg6.mat", ...
           PulseFactor=10, FlowSurge=1.0, Temperature=19, LightScale=1, ...
           TPost=3.0, NFrames=145, NCells=500, Zeta0=1, Zeta1=0.27, DetachForm="linear", ...
           KDOM=3e-4, KHPO4=1e-6, TransferScale=10, EtaSand=1500, Delta=5e-3, MaxDt=5e-5, ...
           Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])
```

### Parameters

Host: **identical to E4** (2× field influent, ζ₀ = 1, ζ₁ = 0.27, linear detachment ×1,
transfer ×10, K_DOM 3e-4, K_HPO₄ 1e-6, η_sand 1500, δ 5 mm, N = 500, MaxDt 5e-5, shin,
Neumann, dark floor 0, no respiration reaction), restarted from
`chain_fld2x_{lit,dark}_leg6.mat` at **t = 60 d** and run 3 further days, 145 frames (30 min).

Pulse: baseline influent PAT **0** (as the snapshots were grown), spiked to
`PulseFactor × 5.36e-3 kg/m³` on t ∈ [60.1, 60.3) d. 5.36e-3 is the Table B.1 marker
influent (`manuscripts/AWR-SSF/results.tex:21`) and the nominal `logOatCampaign.m:35` uses.
`L(t) = log₁₀(C_in,peak / c_PAT,out(t))`.

Four conditions × {lit, dark}: **nom** ×10 pulse; **big** ×100; **surge** ×10 pulse with
inflow velocity ×1.4 for the whole run; **cold** ×10 pulse at **3 °C**.

> **The cold arm is a cold shock, not a winter filter.** Only the θ-corrected kinetics see
> 3 °C; the biofilm it acts on was grown for 60 d at 19 °C. Do not report it as a season.

### Wall time

| task | tag | wall (s) |
|---|---|---|
| `3549855_0` | e4_lit_nom | 601 |
| `3549855_1` | e4_dark_nom | 617 |
| `3549855_2` | e4_lit_big | 737 |
| `3549855_3` | e4_dark_big | 735 |
| `3549855_4` | e4_lit_surge | 646 |
| `3549855_5` | e4_dark_surge | 639 |
| `3549855_6` | e4_lit_cold | 279 |
| `3549855_7` | e4_dark_cold | 286 |

≈ 200 s per simulated day for the 19 °C arms — 2.2× the 92 s/d measured locally, consistent
with the E6/E9 cosmos factor plus seven tasks sharing cn159. The 3 °C arms are 2.2× *faster*
(slower kinetics, larger adaptive steps). The BigPulse arms are ≥ 7× slower than nominal.

### Output and figures

Data: `analysis/probes/data/pulse/pulse_e4_{lit,dark}_{nom,surge,cold}.mat` (133 MB, fetched
by rsync), each with `rec`, `results`, `results_py`. New code:
`analysis/probes/probePulse.m`, `slurm/pulse_e4.sbatch`. **Figures: none.**

### Outcome — removal (unscored) and host state

| tag | L(1 d) | L(3 d) | L_min | peak c_out (kg/m³) | O₂ eff mean / min | O₂ consumed | biomass 0 → 3 d | sup. |
|---|---|---|---|---|---|---|---|---|
| e4_lit_nom | 4.17 | 10.73 | 0.11 | 0.0417 | 5.73 / 5.50 | **3.37** | 0.00797 → 0.00905 (+13.5 %) | 0.000 % |
| e4_dark_nom | 4.18 | 10.74 | 0.10 | 0.0421 | 5.73 / 5.59 | **3.37** | 0.00766 → 0.00867 (+13.2 %) | 0.000 % |
| e4_lit_surge | 4.33 | 10.95 | 0.08 | 0.0448 | 6.39 / 6.07 | 2.71 | 0.00797 → 0.01014 (+27.2 %) | 0.000 % |
| e4_dark_surge | 4.34 | 10.96 | 0.07 | 0.0451 | 6.41 / 6.14 | 2.69 | 0.00766 → 0.00980 (+27.9 %) | 0.000 % |
| e4_lit_cold | 2.83 | 4.82 | 0.11 | 0.0417 | 7.49 / 6.07 | 1.61 | 0.00797 → 0.01158 (+45.3 %) | 0.000 % |
| e4_dark_cold | 2.84 | 4.83 | 0.11 | 0.0422 | 7.48 / 6.14 | 1.62 | 0.00766 → 0.01120 (+46.2 %) | 0.000 % |
| **e4_lit_big** | 3.93 | 10.49 | **0.16** | **0.367** | **1.06 / 0.09** | **8.04** | 0.00797 → 0.02437 (+206 %) | **5.90 %** |
| **e4_dark_big** | 3.94 | 10.50 | **0.16** | **0.372** | **1.05 / 0.17** | **8.05** | 0.00766 → 0.02321 (+203 %) | **3.49 %** |

Units: O₂ mg/L; biomass = ε-weighted ∫φ_b over the whole column, m; "sup." = supernatant
share of it.

**The filter does not remove the pulse while it is arriving.** L_min is **0.07–0.11** at the
breakthrough peak, i.e. peak effluent is 78–84 % of peak influent (0.0417–0.0451 against
0.0536). Removal is entirely a *tail* property: L reaches 4.2 at 1 d and 10.7 at 3 d once the
feed stops. Against the manuscript's claimed ≈ 1.48 log at the end of the feeding time
(`pathogen.tex:76`) this is **an order of magnitude worse at the peak**, on a 60 d mature bed
rather than the manuscript's 3 d one. That gap is the single most important number here and
it is not explained.

**Light does nothing to removal.** lit and dark differ by ≤ 0.01 log at every time. Expected:
the E4 lit/dark bed contrast is 4 %.

**The hydraulic surge makes removal *better*, not worse** (L(3 d) 10.95 vs 10.73, L_min 0.08
vs 0.11 — worse at the peak, better in the tail). ×1.4 velocity raises detachment
(linear law) and cuts residence time, but it also raises the flowing→matrix attachment flux.
Not a robustness result: the surge runs the whole 3 d, not just the pulse.

**Cold halves the removal** (L(3 d) 4.82 vs 10.73). θ-corrected bacterivory at 3 °C is
8.0 × 1.08^(3/293−1) ≈ 3.1 /d against 7.4 /d at 19 °C, and the tail is bacterivory-driven.

**The BigPulse breaks the filter, and that is the most informative run of the eight.** A ×100
spike (peak 0.536 kg/m³ for 4.8 h) drives effluent O₂ from 6.07 to **0.09 mg/L** — consumed
O₂ **8.04 mg/L**, 88 % of the influent, far outside Elemo2024's 2–5 window and below any
aerobic-filter observation. Total biomass **triples** (+206 %), the supernatant fills to
**5.9 %** of it (0.000 % in every other arm), and φ_b at the sand surface reaches **0.67**
(0.37 nominal). Removal at the peak is *worse* than nominal (L_min 0.16 vs 0.11 — the
**relative** breakthrough 0.367/0.536 = 68 % is actually lower than nominal's 78 %, so the
bed is not saturating; the absolute load is simply 10× larger). The lit/dark phototroph
contrast above the sand jumps to **4.38×** — but that is trapped marker plus trapped
influent, not photosynthesis, exactly the artefact `.claude/CRITIQUE.md` §10 identified for
the ζ₀ = 5 arms. **This arm should not be used as a physical prediction**; it is a
demonstration that the marker is dosed as a mass-bearing particle and at ×100 it becomes a
particulate load comparable with the influent solids.

**Every arm gains biomass over the 3 d**, +13 % (nom) to +46 % (cold), against a host that was
flat to +0.1 % per 10 d at 60 d (E9). This is **not** a physical growth spurt — it is the PAT
pulse itself: PAT is a particle and enters φ_b. The cold arms gain most because the marker is
removed slowest. **Consequence: the pulse perturbs the metric it is scored on**, and total
biomass is not a clean host-state observable for a pulse run. Flagged, not fixed.

### Scores (rubric of 2026-08-27; applied by hand — see caveat 1)

Profile is **capped at 3** for every run: the rubric's fourth point needs a 10 d steady window
and these runs are 3 d.

| tag | profile | biomass | O₂ | **total** |
|---|---|---|---|---|
| e4_lit_nom | 3/3 | 3/3 | 3/3 | **9/10** |
| e4_dark_nom | 3/3 | 2/3 | 3/3 | **8/10** |
| e4_lit_surge | 3/3 | 3/3 | 3/3 | **9/10** |
| e4_dark_surge | 3/3 | 2/3 | 3/3 | **8/10** |
| e4_lit_cold | 3/3 | 3/3 | 2/3 | **8/10** |
| e4_dark_cold | 3/3 | 2/3 | 2/3 | **7/10** |
| **e4_lit_big** | **2/3** | **2/3** | **0/3** | **4/10** |
| **e4_dark_big** | **2/3** | **1/3** | **0/3** | **3/10** |

The two BigPulse arms score separately and are itemised after the common justification.

Justification, one line per axis (identical across the six non-BigPulse arms except where noted):

- **Profile 3/3.** φ_b maximum is at **z = 0.0000 m**, the sand surface, in all six (+1).
  Decay below it is monotone apart from a **single** upward step, largest excursion
  2.5e-4 (lit_nom) to 1.06e-2 (dark_cold) of the peak, with one secondary maximum below —
  awarded (+1) as a sub-1 % numerical wiggle at the ε ramp, **not** upward migration; the
  dark_cold arm at 1.06 % is the one borderline case. Supernatant biomass is
  **0.000 %** of the total in all six (+1). Fourth point unavailable at 3 d.
- **Biomass 3/3 lit, 2/3 dark.** Whole-column biomass is monotone in the surge and cold arms
  and dips by 0.07 % once in the two nominal arms before rising — awarded (+1). Top-2 cm
  **46.3 → 51.0 µg C/g wet (lit) / 43.9 → 48.3 (dark)**, i.e. **185 → 204 dry (lit) /
  176 → 193 (dark)**, against Campos2002's 60 at d97: **within a factor 2 on the wet
  convention (+1), out by 3.1–3.4× on the dry one (would be 0)**. f_dry is unresolved
  (E1 caveat 1); both readings are given and the point is awarded on the wet reading, as in
  E9. Convention: µg C/g dry sand = (m/0.02)·0.374/((1−0.4)·2650)·1e6, wet = dry/4 — the
  `scoreRun` convention, which reproduces E4's "45 wet". Phototrophs above the sand,
  lit/dark **2.54× at t = 0 and 2.66× at 3 d** ≥ 1.5× (+1 to the lit arms only; a dark arm
  cannot earn its own contrast point, as in E9).
- **O₂ 3/3 (nom, surge), 2/3 (cold).** Consumed O₂ (influent 9.10 − effluent):
  **3.37 mg/L** nominal and **2.69–2.71** surge, both inside Elemo2024's 2–5 window (+2);
  **1.61–1.62** cold, which is in the 1–2 band (+1). Effluent minimum is **5.50–6.14 mg/L**
  in every arm, ≥ 3 (+1).

**BigPulse justification.**
- **Profile 2/3.** φ_b maximum still at **z = 0** (+1). Decay below it has one upward step of
  **2.7–3.0 %** of the peak — over the 1 % line used above, but still a single step with no
  second maximum higher than the surface, so (+1) is awarded with the excursion stated.
  Supernatant biomass is **5.90 % (lit) / 3.49 % (dark)** of the total, **over the 1 % line
  → 0**. This is the only axis point any arm lost on the profile.
- **Biomass 2/3 lit, 1/3 dark.** Whole-column biomass is **not** monotone — a 1.4 % dip —
  → **0**. Top-2 cm **97.5 µg C/g wet / 390 dry (lit)** and **95.4 / 381 (dark)** against
  Campos's 60: within a factor 2 on the wet convention (+1), 6.4× out on the dry one.
  Phototroph contrast above the sand **4.38×** ≥ 1.5 (+1, lit only) — but see the warning
  above: at ×100 this is trapping, not growth, so the point is arguably unearned.
- **O₂ 0/3.** Consumed **8.04–8.05 mg/L**, outside both the 2–5 window and the 5–6 partial
  band (**0**). Effluent minimum **0.09 mg/L (lit) / 0.17 (dark)**, far below 3 (**0**).

### Caveats

1. **`analysis/scoreRun.m` could not be used.** It expects `probeChain`'s `rec` schema and
   errors with "Subscripted assignment between dissimilar structures" on `probePulse`'s.
   The scores above were computed by hand from `rec` with a scratch script, using
   `scoreRun`'s documented conventions (0.374 kg C/kg, wet = dry/4, ε-weighted integrals),
   verified against E4's published "45 µg C/g wet". They are therefore consistent with
   E9's numbers but were **not** produced by the same code. Either `probePulse` should emit
   a probeChain-shaped `rec` or `scoreRun` should accept both.
2. Total biomass is contaminated by the pulse mass (see above) — the +13…46 % (and +206 % at
   ×100) is marker, not growth. The φ_b-shape and O₂ axes are unaffected at ×10; at ×100 the
   O₂ axis is **not** unaffected, because the trapped marker is itself oxidised.
4. **The tail log removals are not defensible and must not be plotted.** Two independent
   reasons. (i) L(3 d) ≈ 10.7 is a removal of a marker four orders of magnitude below any
   detection limit — a statement about the model's exponential tail, not about a filter.
   (ii) `analysis/testPathogen.m` found a ~5e-5-of-supply non-closure in the ε-weighted PAT
   budget that is **specific to the export term** (retained species close to 3e-11); export
   is 88 % of supply here, so it is ~5.7e-5 of the exported mass, and the export integral is
   dominated by the breakthrough peak at c_out ≈ 4e-2 — an absolute error of that size is
   larger than the entire tail, where c_out ≈ 1e-12. **Only L_min (0.07–0.16) is safe to
   report**, because it sits at the peak where the relative error is negligible. Raised by
   Agent B and accepted; it sharpens rather than softens the finding.

   **The bound, computed on `pulse_e4_lit_nom.mat`.** Supply 8.0399e-2 kg/m²; cumulative
   export 6.0101e-2 (74.8 % of supply); measured non-closure 6.271e-5 of supply
   = **R = 5.042e-6 kg/m²**, i.e. 8.39e-5 of the exported mass. `L(t)` is trustworthy only
   while the marker mass *still to leave the column after t* exceeds R by a comfortable
   factor; taking 10×:

   | t − t_snap (d) | c_out (kg/m³) | L(t) | tail mass after t (kg/m²) | tail / R | defensible |
   |---|---|---|---|---|---|
   | 0.30 | 2.59e-2 | 0.33 | 5.570e-2 | 11048 | yes |
   | 0.50 | 1.53e-2 | 0.54 | 2.342e-3 | 465 | yes |
   | **0.625** | — | **2.94** | 5.04e-5 | **10.0** | **last defensible point** |
   | 1.00 | 3.64e-6 | 4.17 | 3.474e-6 | 0.69 | **no** |
   | 1.50 | 8.34e-8 | 5.81 | 7.969e-8 | 0.02 | no |
   | 3.00 | 1.01e-12 | 10.73 | 5.42e-17 | 0.00 | no |

   **`L` is defensible up to ≈ 2.9 log, reached 0.63 d after the pulse starts (day 60.63).
   `L(1 d) = 4.17` already fails the test (tail/R = 0.69), and `L(3 d) = 10.73` fails it by
   sixteen orders of magnitude.** The `L(1 d)` and `L(3 d)` columns in the outcome table
   above must therefore be read as model diagnostics only, not as removal claims. The
   manuscript's own quoted value of 1.48 log at end-of-feed sits *inside* the defensible
   band, so the comparison in the text is unaffected. `fig:pulse-outflow` plots
   concentrations, not log removals, and is unaffected; `fig:pat-filtration` plots `F(t)`
   and inherits this ceiling. Diagnosing the term itself (comparing the time-integrated
   `q·c_out` against the solver's `ε_face·v_avg(end)·c_flowing(end)`, `simulate.m:921-944`)
   is still owed.
5. The P1 pulses use one snapshot leg (60 d). The manuscript's `fig:pulse-outflow` is a
   **day-30** pulse; that reproduction is P2 below.

---

## P2 — `fig:pulse-outflow` reproduction on the working set (day-30 pulse) — COSMOS + local (2026-08-27)

**Logged** 2026-08-27. **Status:** both panels produced. **Ran on LUNARC cosmos**, job array
**3549870** (`slurm/pulse_fig.sbatch`, array 0–3), from `chain_fld2x_lit_leg3.mat` (t = 30 d),
pulse day 30→32, `TPost` to day 35, N = 500 working set.

Reproduces Manriquez2026 `fig:pat-low-inact` (amplitude 1× vs 100×) and `fig:pat-high-inact`
(the 1× pulse at three rate sets: Experiment 1 d_PAT 0.02 / p_PAT 8.0, Experiment 2 d_PAT 2.0,
Experiment 3 d_PAT 2e-6 / p_PAT 2.0 — the two extra sets are given in the `pathogen.tex`
caption, so no rates were invented). Amplitudes from `manuscriptExperiments.m:366` (Pulse =
c_ref, BigPulse = 100 c_ref).

### Outcome

- **Tasks 0/2/3 (1×, Exp 2, Exp 3) completed** (15–16 min); **task 1 (100×) CLOGGED at day
  30.431**, cell 0 — a 2 d feed at 0.536 kg/m³ delivers ~7.7 kg/m² of marker against a standing
  biomass of ~0.09 kg/m² (~85×) and buries the surface cell 0.43 d into the pulse. The
  cosmos rerun (job 3549875) and a local rerun both CLOGGED at the same time; the truncated
  curve is drawn from `pulse_fig_p100x_ABORTED.mat`.
- **The bed passes the pulse.** Outflow plateaus at ≈ 80 % of the inflow for the whole 2 d feed
  in every rate set (L ≈ 0.1), against the manuscript's ≈ 1.48 log at end of feed. Only the
  **tail** depends on d_PAT/p_PAT: Experiment 3's slow inactivation gives a ~3 d tail, Experiment
  1/2 clear within ~1 d. Same conclusion as P1, on the manuscript's own day-30 protocol.
- **The 100× BigPulse panel is not reproducible on the working set** — the filter clogs. The
  published panel was drawn at ζ₀ = 1e2 with the campaign detachment law on a 3 d filter, a far
  more permissive host.

### Output

Figures: `analysis/results/figures/repro/fig_pulse-outflow_{low-inact,high-inact}.{png,pdf}`
(semilog-y concentrations, dashed rectangular inflow, shaded pulse window, dotted line at the
8.4e-5-of-peak outflow-term resolution from P1 caveat 4). Code: `analysis/plotPulseOutflow.m`.
Data: `analysis/probes/data/pulse/pulse_fig_{p1x,exp2,exp3}.mat` + `pulse_fig_p100x_ABORTED.mat`.
`probePulse` gained a partial-save-on-abort path (as `probeChain` has) and a guard so the
zero-padded post-abort frames no longer break `pulseMetrics`.

---

## E10 — Recreation support runs: winter 90 d, covered 1 %, 20 d scrapes, day-37 pulses — COSMOS (2026-08-28)

**Logged** 2026-08-28. **Status:** all complete. Support runs for the Manriquez2026 figure
recreation (`analysis/recreateFigsManriquez2026.m`, `analysis/results/figures/recreation/`).

| job | tasks | what | wall |
|---|---|---|---|
| 3550889_0 | 1 | `fld2x_winter`: E4 set at **3 °C** with the **manuscript winter light** (probeChain gained `Temperature` and `LightForm` options, smoke-tested at N=100), Days=90 | 2:03:12 |
| 3550889_1 | 1 | `fld2x_cov01`: E4 set at **LightScale=0.01** (the published "covered"), Days=30 | 1:18:50 |
| 3550884 | 4 | scrapes 0/4/8/12 cm rerun with **Days=20** (`repro_scrape.sbatch` gained a `DAYS` env) | ~1:04 each |
| 3550892 | 4 | day-37 pulses from the 30 d snapshot: 1×, **1e-3×** (the published "Pulse 2" — the low arm, not 100×), d_PAT=2, d_PAT=2e-6+p_PAT=2 | ~22 min each |

### Outcome — the seasonal ordering INVERTS on the working set

`scoreRun` at 90 d (`chain_fld2x_{winter,lit}` legs 1–9):

| | winter (3 °C) | summer (19 °C) | ratio |
|---|---|---|---|
| bed ∫εφ_b (m) | **0.02431** | 0.00798 | **3.0×** |
| total column biomass (kg/m²) | **0.2716** | 0.0891 | 3.0× |
| H/H₀ | **1.351** | 1.095 | — |
| O₂ consumed (mg/L) | 2.94 | 2.96 | ≈1 |
| top-2 cm (µg C/g wet) | 46.1 | 41.7 | 1.11 |

**Winter holds 3× the biomass of summer at 90 d and the highest head loss of any run in the
campaign (H/H₀ 1.35), at identical O₂ consumption** — the opposite of the published
fig:seasons-results (summer ≫ winter) and consistent with Bae2023 Table 2's standing stock
(winter ≥ summer) and with the probe-chain seasonal result. Mechanism (hypothesis, not yet
isolated): θ-corrected loss terms (d_HET θ 1.0725, d_PHO θ 1.08, hydrolysis θ 1.0725) fall
faster with cold than growth does, so standing stock accumulates; activity (O₂) stays flat
because it is substrate-limited, not biomass-limited. The s90 study (2026-08-22) separated the
same inversion into kinetics/influent confounds at N=100; this is the first N=500 working-set
pair. **The published summer ≫ winter ordering came from the pre-audit kinetics.**

Also: the covered-1 % arm at 30 d sits between fully dark and lit (surface peak 0.27 vs 0.30
lit), contrast confined to the roughness layer; the 20 d scrape regrowth has **not** converged
back to the unscraped profile, against the published convergence claim; the day-37 pulses
confirm P2's plateau-at-inflow result at both amplitudes (linear in dose; tails per rate set).

### Caveats

1. The winter arm uses the working set otherwise unchanged (2× field influent all year); a real
   seasonal comparison would also vary the influent. Same caveat as the published figure.
2. f_dry unresolved — biomass ratios are convention-free, but the µg C/g values are not.
3. The winter/summer pair differ in BOTH temperature and light curve, as published; the s90
   variants isolate them at N=100 only.

## E11 — Working-set OAT campaign: pulse, flowstep, startup (29 parameters) — COSMOS (2026-09-02)

**Logged** 2026-09-02. **Status:** complete, all three scenarios. The manuscript sensitivity
section's log-OAT campaign (`logOatCampaign.m`) on the audited fixed model, run for all three
disturbance scenarios that Sept-14 needs: pulse, flowstep and startup.

### Commands

```bash
# cosmos, from the repo root (MPC-SSF-manuscript)
sbatch slurm/oat_pulse.sbatch      # job 3562754, array 0-28
sbatch slurm/oat_flowstep.sbatch   # job 3563342, array 0-28
sbatch slurm/oat_startup.sbatch    # job 3563343, array 0-28
```

Each array task runs one parameter's full OAT triple (baseline + ×2/×½ arms, or the constrained
β scheme below) in one `logOatCampaign` call, one task per parameter (not per sign — the
baseline is shared, see the sbatch header rationale):

```matlab
logOatCampaign("pulse", Preset="workingset", Snapshot="chain_fld2x_lit_leg6.mat", ...
    NCells=500, TPost=3.0, MaxDt=5e-5, NFrames=145, ParamFilter="<param>", OutTag="_<param>")
% flowstep: same call, scenario "flowstep"
% startup: scenario "startup", no Snapshot (clean IC; ripening is the disturbance)
```

All three jobs COMPLETED with exit status 0:0 on every one of the 29 array tasks (87 arms
total).

### Parameters

Audited `workingset` preset, anchored on the E4/E9 mature-filter snapshot
`analysis/probes/data/chain/chain_fld2x_lit_leg6.mat` (pulse and flowstep only — startup starts
from a clean column). `NCells=500`, `TPost=3.0` d post-disturbance, `MaxDt=5e-5`, `NFrames=145`.
The 29-parameter design is `campaignParams` inside `analysis/logOatCampaign.m` (velocity
excluded as a fixed operating condition, decision 2026-08-18): `temperature, influent_PAT,
dispersivity, transport_P, attach_sand, sand_pathogen, mu_HET, mu_PHO, d_HET, d_PHO, hydrolysis,
theta_growth, theta_death, K_O2_HET, K_DOM_HET, K_HPO4_HET, K_O2_PAT, K_pred, marker_growth,
inactivation, bacterivory, zeta_0, kappa, zeta_1, detach_scale, light_att_water, light_att_sand,
attenuation_P, beta_porosity` — every parameter perturbed ×2/×½ except `beta_porosity`, which
runs the constrained gap scheme of `analysis/logOatCampaign.m:382-386` (β = 0.98 / 0.95, secant
denominator ln 2.5).

### Wall time

Not pulled from `sacct`; bracketed from the per-parameter output file timestamps (29 tasks run
in parallel per array, 4 h limit each). Pulse: first `L0.csv` 2026-09-01 17:38, last
`measures.csv` 18:36 (~1 h span); per-task time ~31 min (`.claude/JOURNAL.md` 2026-09-01 entry).
Flowstep: first `L0.csv` 23:44, last `measures.csv` 2026-09-02 00:40 (~1 h span). Startup: first
`L0.csv` 23:41, last `measures.csv` 00:15 (~34 min span).

### Output and figures

Data: `analysis/results/oat/log_oat_{pulse,flowstep,startup}_<param>/` (`L0.csv`,
`curves_<param>.csv`, `measures.csv`), one directory per parameter per scenario, 87 total.
Concatenated raw measures: `analysis/results/log_oat_{pulse,flowstep,startup}/measures.csv` (the
pre-existing merged-table convention). Resolution-masked rankings (`maskOatMeasures.m`, LCap =
2.9 log, the O(Δt) closure bound of `2026-09-01-pat-export-closure-resolved.md`):
`analysis/results/oat/measures_masked.csv` (pulse), `measures_masked_flowstep.csv`,
`measures_masked_startup.csv`. **Accepted artefacts for the manuscript** (regenerated
2026-09-02): `analysis/results/figures/oat/oat_ranking_rows_masked.tex` (from
`analysis/results/oat/masked/measures.csv`, the pulse scenario) and
`fig_oat_tornado_masked.{pdf,png}`. The unmasked companions `oat_ranking_rows.tex` and
`fig_oat_tornado.{pdf,png}` are also regenerated but are NOT accepted artefacts — see caveat 1.

### Outcome

**All 87 arms flag OK; `clog_driver = 0` in every one of the three campaigns' 29-row measures
tables** (`log_oat_{pulse,flowstep,startup}_*/measures.csv`) — no OAT arm, at either ×2/×½ or
the constrained β gap, drove the working set to clogging.

Masked ranking leader in all three scenarios is `sand_pathogen` (the PAT-sand sticking
efficiency, block transport), by a wide margin over the rest of the table:

| I_rms (masked) | pulse | flowstep | startup |
|---|---|---|---|
| sand_pathogen | **0.488** | **0.351** | **0.499** |
| attach_sand | 0.157 | 0.097 | 0.021 |
| dispersivity | 0.049 | 0.021 | 0.040 |
| beta_porosity | 0.039 | 0.027 | 0.006 |
| temperature | 0.041 | 0.027 | 0.003 |

(`analysis/results/oat/measures_masked{,_flowstep,_startup}.csv`, sorted by `I_rms_masked`,
top 5 rows of each.) `K_O2_HET` and `K_O2_PAT` are flat zero in all three (audited HET/PAT
Monod terms not limiting on this host, per the 2026-08-27 audit).

### Caveats

1. **The unmasked `oat_ranking_rows.tex`/`fig_oat_tornado.{pdf,png}` are not a pulse/flowstep/
   startup blend — they are the flowstep scenario only.** `readOatMeasures` (shared by
   `makeOatRankingTable.m` and `plotOatTornado.m`) globs `analysis/results/oat/log_oat_*/
   measures.csv`, which matches all three scenarios' per-parameter directories under the same
   parent, then keeps one row per `param` via `unique(T.param, "stable")` — the first hit in
   `dir`'s alphabetical order, i.e. `log_oat_flowstep_*` before `log_oat_pulse_*` before
   `log_oat_startup_*`. Confirmed numerically: `oat_ranking_rows.tex`'s β row (I_rms 0.028) and
   dispersivity row (1.2) match `log_oat_flowstep_beta_porosity`/`_dispersivity` exactly, not
   pulse or startup. The **masked** table does not have this problem — it reads the single
   pre-scoped file `analysis/results/oat/masked/measures.csv`, which is pulse-only by
   construction. Only the masked artefacts should be cited as "the" campaign ranking; the
   unmasked pair is mislabeled and should be regenerated per-scenario before any manuscript use.
2. The masked tables' `kept_fraction` is ~0.94-0.97 (27-28 of 29 disturbance-window samples
   inside the 2.9-log resolution band) — the resolution mask removes only a small tail, not a
   large share of each curve.
3. Startup has no snapshot dependency (clean-column IC), so it is not directly comparable to
   pulse/flowstep's shared `chain_fld2x_lit_leg6` anchor; its ranking is reported alongside the
   other two as a robustness check, not folded into one number.

## Staged (not yet accepted)

| job | script | what | submitted | ETA |
|---|---|---|---|---|
| 3563617 | `slurm/patscrape.sbatch` | fig:pat-scraping re-run — scrape leg3 (30 d) at 0/4/8/12 cm, 20 d regrowth, then constant feed + 10× pulse + 100× pulse via probePulse | 2026-09-02T10:21:16+02:00 | ~2.5 h/arm, 4 arms parallel (supersedes cancelled 3563613) |
