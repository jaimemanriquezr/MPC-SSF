# Development journal

Chronological log, newest entry appended at the bottom. One section per session:
what was done, what was learned, what's next. See `CLAUDE.md` for the convention
and how this relates to `.claude/decisions/`.

---

## 2026-07-28 — Record-keeping setup

### Done

- Established the plans / decisions / journal conventions and documented them in
  `CLAUDE.md` under "Plans and decisions". The flow is mandatory for non-trivial
  work: plan before implementing, decision file when an alternative was genuinely
  considered, journal entry before ending a session.
- Created `.claude/plans/` and `.claude/decisions/` (both empty; the first real
  entries land with the next piece of work).
- Seeded this journal.

Earlier in the session: a status assessment of the project. No code changed.

### Learned

- **`.claude/` was git-ignored in this repo, and is now partially un-ignored.**
  `.gitignore:1` is `.*/`, which matches every dot-directory, so the record-keeping
  paths started out untracked — durable on disk, but unversioned, absent from a
  fresh clone, and invisible to collaborators, which defeats the purpose of
  `decisions/` as the durable searchable "why". Fixed at the end of `.gitignore`.
  Git cannot re-include a path whose parent directory is excluded, so a bare
  `!.claude/plans/` would not have worked; the block un-excludes `.claude/`,
  re-excludes its contents with `.claude/*`, then allows back exactly
  `plans/`, `decisions/`, and `JOURNAL.md`. Verified with `git check-ignore`: the
  three targets and files inside them are trackable, while `settings.local.json`,
  `agents/`, `references/`, `old-src-notes.md`, `julia-port-outline.md`, and
  `references.md` all stay ignored, as do unrelated dot-directories like
  `.vscode/`. Note that the two new directories are empty, and git does not track
  empty directories — they will enter version control with their first file.
- **`CLAUDE.md` was created at the repo root** (`MPC-SSF/CLAUDE.md`), not at the
  `SSF/` parent. The parent file is a stub, sits outside any git repo, and spans
  three sibling repos, so `.claude/plans/` would have been ambiguous there. At the
  repo root the paths resolve unambiguously and the file is trackable.

### Next

- **The Manriquez2026 revision is the live deadline: Aug 1, 2026.** The sensitivity
  analysis the reviewers asked for is computationally complete (Log-OAT across three
  scenarios, Sobol, mesh-convergence and mesh-robustness checks, reference data from
  Campos2006 / Schijven2013) and documented in `julia/analysis/SENSITIVITY_ANALYSIS.md`
  and `SESSION_LOG_2026-07.md`. What does not exist is the deliverable: no
  response-to-reviewers text, no manuscript subsection, no formatted figure. The
  remaining work is writing, not computing.
- Coverage gap against the reviewer list in `TODO.md`: light attenuation
  (`light_attenuation_{water,sand}`, particle `attenuation`) and the cohesion
  parameters `kappa` / `zeta_1` are not among the 22 Log-OAT parameters. Light was
  swept in the earlier `julia/analysis/sensitivity.jl` run and came out near-zero,
  with a physical explanation (at `light_attenuation_sand=1500` light penetrates
  0 mm into the sand, so the effect is structurally absent) — but that argument is
  not currently backed by a number in a results file. R1 asked about it specifically.
- `TODO.md`'s "Article revision" section still has all five boxes unchecked though
  four are effectively done; it is misleading as a status source.

---

## 2026-07-28 (second session) — light attenuation added to the Log-OAT sweep

### Done

- **Closed the light-attenuation coverage gap** flagged as "Next" in the previous entry. Added
  three parameters to `analysis/log_oat_sensitivity.jl`, taking `PARAMS` from 22 to 25:
  `light_att_water` (0.32, `SandFilter.light_attenuation_water`), `light_att_sand` (1500.0,
  `SandFilter.light_attenuation_sand`), and `attenuation_P` (0.094, `Particle.attenuation`
  across all four particles). New `"light"` block.
- Setters follow the two patterns already established in the file: `f`-mutating for filter
  fields (as `set_velocity`/`set_temp` do — safe because `f2` is freshly constructed per run at
  `log_oat_sensitivity.jl:159`), and `_mdl`/`_pt` rebuilding for the particle field (as
  `set_disp` does).
- Ran `main(disturbance=:startup)` — 51 runs. Log at
  `../../logs/log_oat_startup_2026-07-28.log`; outputs overwrite
  `analysis/results/log_oat_startup/` (the previous 22-row version is recoverable from git —
  114 result files are tracked and were clean before the run).
- Plan: `.claude/plans/2026-07-28-light-attenuation-oat.md`. Decision:
  `.claude/decisions/2026-07-28-light-oat-startup-scenario.md`. These are the first real entries
  in both directories, so they also validate that the un-ignore block added to `.gitignore` last
  session actually makes them trackable.
- Wrote a co-author brief at `../../reports/sensitivity-analysis-coauthor-brief.typ` (compiles
  clean). Jaime's near-term deliverable is a co-author meeting, not the rebuttal.

### Learned

- **The `:pulse` scenario is structurally incapable of measuring light sensitivity.**
  `log_oat_sensitivity.jl:190` builds the mature state **once** under the nominal model, and every
  perturbed run re-homes that frozen state (`:171-172`) while the perturbation touches only
  `(f2, m2)` (`:166`) — i.e. only the 1.5-day challenge window. Light acts *cumulatively* through
  phototroph growth, so perturbing it after the biomass is fixed reports `I_rms ≈ 0` as an
  **artifact of the protocol**. That matters precisely because near-zero is the answer we want to
  give R1: it is only defensible from a scenario where non-zero was reachable. `:startup` uses a
  clean IC (`:170`), so ripening happens inside the run. This generalizes — **any parameter whose
  effect is cumulative rather than instantaneous is mis-measured by the `:pulse`/`:flowstep`
  protocol**, which is worth checking before adding further parameters.
- The two filter-side light knobs need separate treatment, not lumping: `light_att_sand=1500`
  extinguishes light within a fraction of a mm of the bed surface and still does at 750, so it is
  inert **by construction** (`src/SandFilter.jl:112-123`) — an argument independent of any
  measurement. `light_att_water=0.32` governs the supernatant
  (`src/SandFilter.jl:106-109`), where the manuscript's own abstract says biofilm grows "up into
  the supernatant water", so it is the informative case.
- **Presentation hazard in the existing results**: `velocity` has the largest `I_rms` (0.763) in
  `:pulse` but is flagged `clog_driver` — it clogs at ×2, as do `beta_porosity` (0.337) and
  `zeta_0` (0.071) at ×½. Their elasticities are regime changes, not smooth removal
  sensitivities. Reporting one ranking over both kinds would overstate `velocity` as the most
  influential parameter when what it does at ×2 is terminate the run.
- **The Sobol run is weaker than the table suggests**: at N=64 the first-order `S_i` CIs straddle
  zero for both top parameters (`sand_pathogen` `S_i=0.279`, CI −0.113 to 0.904). Only `S_T` is
  resolved, and `S_T` summing above 1 for the pair indicates real interaction between the two
  attachment terms. Report `S_T` only.
- **MATLAB side cannot host this analysis.** No pathogen model exists in `src/` on any branch
  (`main`, `dev`, `matlab-claude`, `julia-port`) — the only hit is a comment at
  `matlab-claude:src/@State/simulate.m:290`, and the only preset is `modelLund`. Also missing:
  `implicit_osmosis` (the ~20× speedup), mature-state resume, and any sweep driver. Separately,
  `Reaction.EfficiencyFlowing` **is declared but never read** by `@State/simulate.m` or
  `@Model/computeReactionRates.m` — the `water_factor` analogue is not wired into the solver.
  `TODO.md`'s unchecked adaptive-CFL box is **stale for `matlab-claude`**, which does have it.

### Next

- Insert the measured light elasticities into the co-author brief and regenerate figures with
  `analysis/render_log_oat.jl` against `results/log_oat_startup/`.
- Re-run `:pulse` and `:flowstep` at 25 parameters for a complete cross-scenario table.
- **Run `Pkg.test`** — not done since `log_oat_sensitivity.jl` was edited. The analysis driver is
  not part of the golden suites, so this is expected-safe but unverified.
- Fix `TODO.md`: tick the four effectively-done "Article revision" boxes, and correct the stale
  adaptive-CFL box for `matlab-claude`.

---

## 2026-07-28 (third session) — cohesion potential fixed; MATLAB pathogen port; cosmos set up

### Done

- **Fixed the Cahn-Hilliard potential to the published form.** `CahnHilliardModel`'s default
  `potential_gradient` was `0.25·u²(1−u)²`; the manuscript and **all three** legacy `.mat` models
  implement `Ψ′ = u²(u − 3ζ₁/2)`. Changed the default and made `zeta_1` live for the first time.
  Reasoning and the full evidence table are in
  `.claude/decisions/2026-07-28-published-cohesion-potential.md`; plan in
  `.claude/plans/2026-07-28-fix-cohesion-potential.md`.
- Added `set_kappa()` and `set_zeta1()` to `analysis/log_oat_sensitivity.jl` and registered `kappa`
  (1e-6) and `zeta_1` (1e-2) in `PARAMS`, now 27 parameters. Rewrote `set_zeta0()` to carry the
  cohesion handles through a new `_ch` helper.
- Updated the `CahnHilliardModel` unit test, which asserted the old form, and extended it to check
  the legacy handle values and that ζ₁ actually moves the potential.
- **cosmos is set up and working.** `~/MPC-SSF` populated by rsync; reusable sbatch template at
  `slurm/pkgtest.sbatch`; logs under `slurm/logs/`. Verified: Julia 1.10.4 instantiates the project
  (stdlib-only deps, `compat julia = "1.10"`), and both Symbolic Math and Parallel Computing
  toolboxes are licensed — the former is effectively required by the legacy `simulate_filter`.
- **MATLAB pathogen port landed and committed** on `matlab-claude` as **`b37e533`** (started by a
  worktree agent, finished here after it hit a session limit): added
  `Particle.SandAttachmentFactor`; wired the previously-dead `Reaction.EfficiencyFlowing` /
  `EfficiencyBiofilm` into the three `evaluateReactions` calls; split the flowing attachment term
  into sand/biofilm form; added `src/presets/modelPathogen.m` and
  `julia/test/golden/compare_pathogen.m`.
- **MATLAB cohesion potential fixed**, same commit `b37e533`. Found a *second* instance of the
  closure trap: `Model.m` assigns `Zeta1` field-by-field **after** constructing the submodel, which
  left the handle built from the old value, so `Model.m` now rebuilds through the constructor.
  MATLAB cannot reference a sibling argument in an `arguments` default, so `PotentialGradient` is
  declared with **no** default and built in the constructor body, using `isfield` to detect
  omission.
- **Light golden coverage added** as **`4ffe5bb`**. `compare_pathogen.m` had hardcoded
  `DarkRespiration=0` and the uniform dark seed, so `reference_pathogen_light` was unvalidated on
  the MATLAB side; it now takes `DarkRespiration`, `LightConst`, `SeedMatrix` and
  `SeedEnclosedLiquids`, mirroring `PGOLD_LIGHT_RUN` in `compare_pathogen.jl`.
- **Sensitivity sweeps re-running on cosmos** as job array **3427577** (3 tasks: pulse, startup,
  flowstep; 27 parameters, 55 runs each; nodes cn127/cn169/cn018). Reusable template at
  `slurm/logoat.sbatch`.

### Verification

- **216/216 tests pass** on cosmos (job 3427570, Julia 1.10.4, 26 s), including all four
  golden-master suites.
- The new default reproduces the legacy handles to machine precision: worst absolute difference
  **0.0** vs `SDparameters*.mat` (ζ₁=0.005) and **2.2e-19** vs `thesis_model.mat` (ζ₁=0.01).
- Setter behaviour asserted directly: `set_zeta1(0.02)` moves `dpsi(0.05)` to 5.0e-5 while
  `set_zeta0` and `set_kappa` leave it at the baseline 8.75e-5.
- **`modelLund` is bit-identical under the MATLAB port, by construction not by measurement**: the new
  `attachmentFlowingFactorP = (1 − poro).*sandFactors + poro.*phiBiofilm` with `sandFactors = 1` is
  algebraically the old `attachmentFlowingFactor = (1 − poro) + poro.*phiBiofilm` (`simulate.m:288`),
  and `EfficiencyBiofilm`/`EfficiencyFlowing` default to exactly `1.0`. Multiplication by 1.0 is
  exact in IEEE arithmetic. `Particle`'s constructor guards `if ~isempty(input_value)`, so the 1.0
  property default survives when the argument is omitted.

### Learned

- **The goldens could never have caught the potential bug, and still cannot.** `reference/` and
  `reference_adaptive/` start from clean ICs where φ_b = 0, and both potential forms vanish at 0;
  the pathogen goldens seed φ=0.05 but run only 1e-5 d. Confirmed empirically — all four suites pass
  *unchanged* after a fix that alters the potential by 98× at φ=0.01 and 18× at φ=0.63. A
  discriminating golden (seeded biofilm, long horizon) is still owed.
- The wrong form was not merely offset but **sign-changing across the physical range**: too large
  below φ≈0.2, too small above it. At the documented self-limiting peak φ=0.63 cohesive resistance
  was ~18× too weak.
- **Closure trap**: because the default closes over `zeta_1`, passing `potential_gradient` through a
  reconstruction carries the *old* ζ₁. A naively written `set_zeta1` would change the field and leave
  physics untouched — a silent zero sensitivity that looks like a result.
- Two cosmos gotchas: sbatch scripts need `#!/bin/bash -l` or `module` is not found; and nested
  quoting through `ssh` + `matlab -batch` is unusable — write the MATLAB to a file via heredoc.

### Verification (MATLAB side)

- 5/5 cohesion checks pass: degenerate `u^3` when `Zeta1` absent; **zero** absolute difference
  against the `SDparameters*.mat` handle; an explicitly supplied handle still honoured; and
  `modelLund` / `modelPathogen` both give `dpsi(0.05) = 8.75e-05`, matching Julia — which is the
  check that proves the post-construction `Zeta1` assignment now rebuilds the closure.
- Both pathogen goldens match: dark at `atol=1e-7` (worst **6.37e-09**), light at `atol=1e-6`
  (worst **2.66e-07**). The relaxed tolerance for the light case is the same one the Julia side
  needs (`test/golden/README.md` records worst ~3e-07), and the only two fields above 1e-7 are
  `PHO/matrix` and `O2/enclosed` — phototroph biomass and photosynthetic oxygen, where
  cross-implementation float divergence is expected to concentrate.

### Correction to an earlier claim in this entry's session

I reported mid-session that the light parameters caused a ~10x solver slowdown. **That was a
measurement artifact and is withdrawn.** `log_oat_startup/` and `log_oat_flowstep/` already held 22
committed `curves_*.csv` files, so counting files in those directories showed "22/27" regardless of
progress — the job was overwriting them one at a time. Real progress must be read from the
newest-modified file, or from `log_oat_pulse/`, which is a fresh directory. Whether light is
expensive is **still an open question**. Measured rate across all three scenarios is ~3.7 min per
parameter.

### Next

- **Commit the Julia-side changes.** `julia/src/cohesion/CahnHilliardModel.jl`,
  `julia/analysis/log_oat_sensitivity.jl` and `julia/test/runtests.jl` are **uncommitted on
  `julia-port`**, and the cosmos sweeps are running against them. This is the top risk.
- Add the discriminating cohesion golden — still missing. All four suites passed unchanged across a
  98x change in the potential, so nothing currently tests it.
- When array 3427577 finishes, diff the new 27-parameter rankings against the committed
  `measures.csv` to quantify what the corrected cohesion law moved.
- Biomass-QoI steady-state study (workstream B) not started. The criterion had to be relaxed to
  1%/day: biofilm mass follows a power law t^0.85 and never plateaus, so 0.1%/day is ~850 simulated
  days away.
- `matlab-claude` has **no upstream**, so `b37e533` and `4ffe5bb` are local only.
