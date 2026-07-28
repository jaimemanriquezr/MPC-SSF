# Plan: add light-attenuation parameters to the Log-OAT sweep

## Goal

Reviewer R1 of Manriquez2026 asked specifically about light attenuation. The 22 parameters in
`analysis/log_oat_sensitivity.jl:PARAMS` do not include any light knob, so the revision
currently has no number to cite — only a verbal argument that the effect is negligible.

**Done** means: the three light parameters are swept by the existing Log-OAT machinery, their
elasticities are written to a durable results file under `analysis/results/`, and the number can
be quoted in the co-author report. Confirming the effect is negligible is a *fine* outcome; the
point is that it becomes measured rather than asserted.

The three knobs (all present, all with nominal values):

| Parameter | Nominal | Where |
|---|---|---|
| `light_att_water` | 0.32 | `SandFilter.light_attenuation_water` (supernatant optical depth) |
| `light_att_sand` | 1500.0 | `SandFilter.light_attenuation_sand` (roughness layer + packed bed) |
| `attenuation_P` | 0.094 | `Particle.attenuation`, all four particles in `modelPathogen` (biofilm self-shading) |

## Protocol subtlety — this drives the design

Light acts on the phototroph population **cumulatively**, through growth. But
`log_oat_sensitivity.jl:190` builds the mature state **once, under the nominal model**, and every
perturbed run re-homes that same frozen state (lines 171–172). A light perturbation applied only
to the 1.5-day challenge window therefore *cannot* move the phototroph biomass that light
governs — it would report ~0 **because of the protocol, not because of the physics**.

Reporting that as "light doesn't matter" would be a false negative, and it is exactly the kind of
claim a reviewer can puncture.

The `:startup` disturbance is the honest home for this question: it uses a clean IC
(`log_oat_sensitivity.jl:170`), so ripening happens *inside* the run and light is free to act
across the whole simulation.

**Decision**: sweep the light parameters in **all three** scenarios for consistency with the
other 22, but treat `:startup` as the scenario the reviewer response cites, and state the
`:pulse`/`:flowstep` limitation explicitly rather than burying it. Recorded in
`.claude/decisions/2026-07-28-light-oat-startup-scenario.md`.

Physical expectation to test (not assume): `light_att_sand=1500` gives an optical depth that
extinguishes light within ~0 mm of the bed surface, and halving it to 750 should still
extinguish — so `light_att_sand` should be structurally inert. `light_att_water` (0.32) governs
the **supernatant**, where the manuscript's own abstract says biofilm grows "up into the
supernatant water" — so this one is *not* obviously negligible and is the interesting case.

## Steps

1. Add three setter closures beside the existing ones (`log_oat_sensitivity.jl:102–116`),
   following the established two patterns exactly — `f`-mutating for filter fields (as
   `set_velocity`/`set_temp` do; `f2` is freshly constructed per run at line 159, so in-place
   mutation is safe), and `_mdl`/`_pt` rebuilding for the particle field (as `set_disp` does).
   - `set_light_water()`, `set_light_sand()` → mutate `f`
   - `set_attenuation()` → `_pt(c; attenuation=v)` over all `Particle`s
   *Verify*: `julia -e` loads the file without error; `length(PARAMS) == 25`.

2. Add three `P(...)` rows to `PARAMS` in a new `"light"` block, with `source` provenance
   consistent with neighbouring entries.
   *Verify*: the three names appear in `measures.csv` with `block == "light"`.

3. Run `main(disturbance=:startup)`, then `:pulse` and `:flowstep`.
   *Verify*: each `results/log_oat_<scenario>/measures.csv` has 25 rows; three new
   `curves_light_*.csv` / `curves_attenuation_P.csv` files exist per scenario.

4. Record the measured `I_rms` for each light parameter and its rank among 25.
   *Verify*: numbers land in a durable file, not only in conversation (per `CLAUDE.md`).

## Files touched

- `julia/analysis/log_oat_sensitivity.jl` — three setters + three `PARAMS` rows. **Note this
  file already carries uncommitted user changes (4 insertions); preserve them.**
- `julia/analysis/results/log_oat_{startup,pulse,flowstep}/` — regenerated outputs.
- `.claude/decisions/2026-07-28-light-oat-startup-scenario.md` — new.

## Risk

Regenerating `measures.csv` overwrites the existing 22-row results. Those are committed
(`git ls-files julia/analysis/` = 156 files), so they are recoverable via git — confirm before
the first run.
