# Wet/dry mass convention: f_dry = 0.25 adopted, as a side-project (2026-09-01)

## The decision

Jaime, 2026-09-01: "Adopt f_dry = 0.25 and submit the OAT on the current anchors,
but as a side-project."

Three parts:

1. **f_dry = 0.25 is the model's dry-mass fraction of particulate biomass** — the
   §G fix (state variable is wet cell mass, RWQM1 stoichiometry is per dry organic
   matter; yields were ~4× too strong per kg of state variable). Fix as specified in
   JOURNAL 2026-08-26: ρ_P/f_dry with influents scaled; φ_b, head loss, clogging and
   detachment invariant — that invariance is the verification gate.
2. **The manuscript OAT runs on the current anchors** (pre-f_dry convention),
   submitted 2026-09-01 as cosmos job 3562754 (`slurm/oat_pulse.sbatch`, array 0-28,
   E4 snapshot `chain_fld2x_lit_leg6`). Deliberate: the OAT's claim is a *relative*
   log-sensitivity ranking, robust to a uniform yield rescale, and re-anchoring
   first would cost days the Sept-14 deadline does not have.
3. **Side-project status**: implementing f_dry in the presets and regenerating the
   anchors (E-chains, 104 d pair, seasons/roofed/scrape data) happens after the
   Friday co-author draft, targeting the Sept 14 resubmission or the response
   letter. Until then the parameter table declares biomass densities as wet-mass
   values.

## Why 0.25

- §G's own "~4×" estimate implies f_dry ≈ 0.25.
- Melo2005 (`documents/Melo2005.pdf`, Wat. Sci. Tech. 52(7):77-84) makes it citable
  at the correct scale: via Hinson & Kocher 1996 (cited p. 80), dry cell material is
  200-670 kg/m³ against ~1000 kg/m³ hydrated cells → cell-scale f_dry ≈ 0.2-0.67;
  0.25 is the credible low end. The *biofilm-scale* 90-99 % water content (f_dry
  0.01-0.10) is NOT the right number here: the model already carries that water
  structurally as the enclosed-water phase (β = 0.99); using it would double-count.

## Consequence to carry: Melo couples f_dry to β

Melo Table 1: measured biofilm dry density ρ_dw = 14-91 kg dry/m³ wet biofilm
(nine labs, porosities 0.92-0.99). Model: ρ_dw = (1−β)·ρ_P·f_dry. With β = 0.99,
ρ_P = 1117, f_dry = 0.25 → 2.8 kg/m³, below every entry (even f_dry = 1 gives 11).
Matching the table's low end needs β ≈ 0.95; typical 30-40 kg/m³ needs β ≈ 0.87-0.91.
**When f_dry is implemented, β = 0.99 should be revisited** — the sensitivity draft
holds β at 0.99 (26→25 parameters), but `oat_pulse.sbatch` still perturbs
`beta_porosity` (parameter 29), so the campaign quantifies the stakes either way.
Caveat: Melo's own conclusion is that ρ_dw correlates poorly with *transport*
(tortuosity dominates) — this constraint is mass bookkeeping only, which is exactly
what the §G yield question is.

## Alternatives rejected

- **Defer f_dry entirely** (recommendation of the 2026-09-01 brief, option 1):
  overtaken by Jaime's decision; the middle path (option 3) is what was adopted.
- **Re-anchor before the OAT**: invalidates every snapshot behind the new figures
  during draft week.
- **Biofilm-scale f_dry (0.01-0.10)**: wrong scale; double-counts enclosed water.

## Files touched (when implemented — not yet)

- `src/presets/modelLund.m`, `src/presets/pathogenModel.m` (ρ_P/f_dry, influent scaling)
- `analysis/PARAMETERS.md` (document the convention)
- verification: φ_b path invariance on a reference run before/after
