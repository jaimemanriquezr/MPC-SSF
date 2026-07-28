# Light-attenuation OAT is cited from the `:startup` scenario, not `:pulse`

## The decision

The three light parameters added to `analysis/log_oat_sensitivity.jl:PARAMS`
(`light_att_water`, `light_att_sand`, `attenuation_P`) are swept in all three disturbance
scenarios for consistency with the other 22, but the **reviewer response cites `:startup`**, and
the `:pulse` / `:flowstep` light numbers are reported with an explicit caveat rather than as
evidence of insensitivity.

## Why

Light acts on the phototroph population **cumulatively**, via growth. The OAT protocol builds the
mature state once, under the nominal model, and re-homes that same frozen state into every
perturbed run:

- `analysis/log_oat_sensitivity.jl:190` — `ms = build_mature(m0, ncells, tmature)`, called once
  with the unperturbed `m0`.
- `analysis/log_oat_sensitivity.jl:171-172` — each perturbed run constructs its `State` from
  `deepcopy(ms.global_concentration)`.
- `analysis/log_oat_sensitivity.jl:166` — the perturbation is applied to `(f2, m2)` only, i.e.
  only to the 1.5-day challenge window.

So in `:pulse` and `:flowstep`, a light perturbation is applied *after* all the biomass that light
governs has already been fixed. Such a run is **structurally incapable** of showing a light
effect. It would report `I_rms ≈ 0` as an artifact of the experimental design.

That distinction matters because the near-zero result is precisely what we want to tell R1, and a
reviewer is entitled to ask whether the protocol could have detected an effect had one existed.
An answer of "we measured zero" is only defensible from a scenario where a non-zero result was
reachable. `:startup` uses a clean initial condition
(`analysis/log_oat_sensitivity.jl:170` — `State(f2, m2)`), so ripening happens inside the run and
light acts across the whole simulation.

Secondary reason to expect a genuine split between the two filter-side knobs, which is why they
are swept separately rather than as one lumped "light" factor:

- `light_att_sand = 1500` enters as `f.light_attenuation_sand .* eta` with
  `eta ∝ (1-ε₀)·(δ/2 + z)` (`src/SandFilter.jl:112-123`). At that magnitude light is extinguished
  within ~0 mm of the bed surface, and halving it to 750 still extinguishes. Expected inert **by
  construction**, and that is a physical argument that stands independent of the number.
- `light_att_water = 0.32` governs the **supernatant** (`src/SandFilter.jl:106-109`,
  `η = 0.32·(z + height)`), where the manuscript's own abstract states biofilm grows "up into the
  supernatant water". This one is not obviously negligible and is the informative case.

## Alternatives rejected

**Rebuild the mature state under each perturbed light value.** This is the most physically
faithful option: it would let `:pulse` report a real cumulative light sensitivity. Rejected for
now on cost and comparability. Cost: two extra 3-day maturation builds per light parameter, and
`build_mature` is the most expensive single call in the harness. Comparability: it would make the
three light rows incommensurable with the other 22, which are all challenge-window perturbations,
so they could not share a ranking table — defeating the purpose of putting them in the same
`measures.csv`. `:startup` gets the same physics for free, because there the maturation *is* the
run.

**Sweep one lumped light factor instead of three knobs.** Rejected: it would hide the fact that
the sand-bed insensitivity is structural (extinction) while any supernatant sensitivity is
physical. Those are different arguments to a reviewer and need different numbers.

**Argue the point in prose with no new computation.** This was on the table and Jaime explicitly
chose to run the sweep and cite a measured number. The prose argument (optical depth extinguishes
light in the bed) remains correct and is worth keeping, but it now leads into a number rather
than standing alone.

## Numbers

Measured elasticities live in `analysis/results/log_oat_startup/measures.csv` (rows with
`block == "light"`), with per-parameter response curves in the sibling `curves_*.csv`. Run log:
`../../logs/log_oat_startup_2026-07-28.log`. This file deliberately does not restate the values —
see the results file, which is the durable record.
