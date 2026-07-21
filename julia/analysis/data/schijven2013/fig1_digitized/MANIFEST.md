# Fig. 1 (Schijven et al., 2013) — digitized measured data points

Automated pixel digitization of the **measured** data points in Figure 1 of
Schijven, Bradford, Yang (2013), *Water Research* 47, 2592–2602 (figure on
**page 6** of `Schijven2013.pdf`).

Figure caption: *"Breakthrough curves: Measured seeding concentrations (filled
small symbols), measured breakthrough concentrations (open large symbols) and
fitted two-site kinetic model using Hydrus-1D (lines) of all experiments with
MS2 and ECWR1."*

Only the **measured markers** were digitized. The fitted Hydrus-1D **lines were
NOT digitized** — they are regenerable from the paper's two-site parameters
already extracted under `../` (Table 2/3).

## Files
One CSV per panel (4×2 grid): rows W/L/D/G × columns MS2 (virus, left) /
ECWR1 (*E. coli*, right):
`W_MS2.csv W_ECWR1.csv L_MS2.csv L_ECWR1.csv D_MS2.csv D_ECWR1.csv G_MS2.csv G_ECWR1.csv`

### Columns
- `panel` – panel id (e.g. `W_MS2`)
- `organism` – `MS2` or `ECWR1`
- `sand_group` – `W`/`L`/`D`/`G`
- `day` – x value (days), from pixel→day calibration
- `concentration` – C, absolute value (not log), 4 significant figures
- `series_type` – `breakthrough` (open large symbol = effluent) or `seeding`
  (filled small gray symbol = influent)
- `experiment` – **best-effort** experiment id (e.g. `W1`), inferred from marker
  *shape* via the panel legend. Blank when shape could not be resolved. **Low
  reliability — see limitations.**
- `confidence` – `high` / `medium` / `low` (position reliability; see below)
- `shape` – detected marker shape (`circle`/`square`/`triangle`/`diamond`/`filled`/`unknown`)
- `units` – `pfp/ml` (MS2) or `cfp/ml` (ECWR1)

## Method
1. **Render**: `pdftoppm -f 6 -l 6 -r 600 -png` → 4961×6615 px grayscale.
2. **Axis calibration** (per panel). The panels are L-shaped (left y-axis +
   bottom x-axis, no top/right frame). The y-axis line (longest vertical dark
   run) gives day-0 / origin x; the x-axis line (longest horizontal dark run)
   gives the bottom-decade y. Major x-tick marks (below the x-axis) and labeled
   y-decade tick marks (left of the y-axis, ticks every 2 decades) were detected
   as evenly-spaced dark segments; `px_per_day` and `px_per_decade` come from
   their median spacing. Bottom-decade *values* were read directly from
   rendered crops of the y-axis labels. Mapping:
   `day = (px_x − x0)/px_per_day`, `log10 C = bot + (y_bottom − px_y)/px_per_decade`.
3. **Open (breakthrough) markers**: detected as **enclosed white holes** — white
   connected components not touching the plot border, of marker-interior size.
   This isolates the open circle/square/triangle/diamond centers and naturally
   ignores the fitted lines (which enclose no white) and the filled gray markers
   (no white center).
4. **Filled (seeding) markers**: detected as solid mid-gray blobs
   (120 < intensity < 208) of single-symbol size.
5. **Shape → experiment**: hole shape classified by fill-extent (area/bbox) and
   top/bottom taper, restricted to the shapes present in each panel's legend.
6. **Legend + title exclusion**: the top-right legend box and top-center title
   band of each panel are masked so legend symbols and letter-holes (e.g. the
   loop in "R") are not mistaken for data.
7. **Confidence** (position): each breakthrough point's neighbor count within
   0.30 day and 0.35 log10-units — 0 neighbors→`high`, 1–2→`medium`, ≥3→`low`.
   This automatically flags the crowded day 0–1.5 peak cloud as low confidence
   and well-separated tail points as high. All `seeding` points are `low`.

DPI: **600**. Libraries: numpy, scipy.ndimage, Pillow (opencv-python-headless
also installed). `pip install --user` succeeded.

## Per-panel calibration constants
(x0 = origin px; ppd = px/day; y_bottom px; ppdec = px/decade; bottom decade value)

| panel | x0 | ppd | y_bottom | ppdec | bottom C | x-tick interval |
|-------|----|-----|----------|-------|----------|-----------------|
| W_MS2   | 1025 | 155.75 | 1491 | 99.5   | 1e-3 | 2 |
| W_ECWR1 | 2876 | 116.0  | 1493 | 99.5\* | 1e-3 | 2 |
| L_MS2   | 1056 | 122.75 | 2618 | 110.25 | 1e-2 (top 1e6) | 2 |
| L_ECWR1 | 2878 | 112.25 | 2620 | 98.0   | 1e-3 | 2 |
| D_MS2   | 1053 | 273.5  | 3753 | 98.0   | 1e-3 | 1 |
| D_ECWR1 | 2877 | 273.0  | 3754 | 98.0   | 1e-3 | 1 |
| G_MS2   | 1053 | 273.5  | 4888 | 98.0   | 1e-3 | 1 |
| G_ECWR1 | 2877 | 339.0  | 4889 | 98.0   | 1e-3 | 1 |

\* W_ECWR1 y-tick detection was noisy (title/label interference); its y-scale
was taken from W_MS2, which is the same figure row (identical y-axis, both
1e-3…1e5).

Note: legends read from the figure — **L-MS2 = L1/L2/L3/L4** but **L-ECWR1 =
L2/L3/L4/L5** (no L1); **D = D1/D2/D3**; **G = G1/G2**; **W = W1/W2/W3**.
L_MS2 uses a taller decade range (0.01…1,000,000) than the other panels
(0.001…100,000).

## Point counts & confidence breakdown

| panel | total | breakthrough | seeding | high | medium | low | experiment-labeled |
|-------|------:|-------------:|--------:|-----:|-------:|----:|-------------------:|
| W_MS2   |  74 |  60 | 14 |  8 | 31 | 21 | 59 |
| W_ECWR1 |  56 |  45 | 11 |  9 | 25 | 11 | 45 |
| L_MS2   | 190 | 161 | 29 | 16 | 94 | 51 | 161 |
| L_ECWR1 | 128 | 116 | 12 | 17 | 66 | 33 | 116 |
| D_MS2   |  79 |  63 | 16 |  3 | 22 | 38 | 60 |
| D_ECWR1 |  94 |  79 | 15 |  3 | 34 | 42 | 78 |
| G_MS2   |  58 |  46 | 12 |  2 | 19 | 25 | 39 |
| G_ECWR1 |  65 |  50 | 15 |  0 | 17 | 33 | 38 |
| **all** | **744** | **620** | **124** | 58 | 308 | 254 | — |

## Accuracy limitations — read before use
This is digitization of a **dense, log-scale, multi-series** figure. Numbers are
approximate; do not treat them as exact.

- **Log y-axis**: a few-pixel error is a multiplicative error in C. Typical
  positional error is roughly ±0.03–0.08 decade (±8–20 % in C) for clean points,
  larger where symbols overlap.
- **Peak cloud (day ≈ 0–1.5)**: all series pile up near the influent
  concentration and symbols heavily overlap. Points here are flagged
  `confidence=low`; counts are **undercounts** (overlapping markers merge into
  one detected blob) and individual (day, C) values are unreliable.
- **Tails (well-separated)**: `confidence=high`/`medium` points, mostly day ≳ 2,
  are the **trustworthy** part of this dataset — the four (L) / three (D) / two
  (G) declining curves are cleanly separated there.
- **Experiment attribution is best-effort and LOW reliability.** It is inferred
  purely from marker shape, which the classifier confuses in crowded regions and
  sometimes even in the tails (circle↔diamond↔square↔triangle). Trust the
  `series_type` (breakthrough vs seeding) far more than the `experiment` label.
  Where you need per-experiment (W1 vs W2 vs W3) curves, prefer well-separated
  tail points and verify shape against the figure. Do **not** assume the
  experiment labels are correct.
- **Seeding (filled) points**: captured where separable; in the densest clouds
  they merge and are undercounted. All marked `low`. Seeding and early
  breakthrough overlap in the same top region.
- **Fitted-model lines were intentionally not digitized.**

### Per-panel trust summary
- **L_MS2, L_ECWR1, D_ECWR1**: tails excellent — 3–4 curves cleanly separated
  and captured; shape/experiment mostly right in the tail. Peak cloud unreliable.
- **W_MS2, W_ECWR1**: good tails (3 curves); moderate peak overlap.
- **D_MS2, G_MS2, G_ECWR1**: fewer high-confidence points (curves plateau close
  together / shorter tails); positions in tails are reliable but experiment
  attribution is weaker. G panels have a large low-confidence fraction.
