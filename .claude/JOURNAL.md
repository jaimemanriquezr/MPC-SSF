# Development journal

Chronological log, newest entry appended at the bottom. One section per session:
what was done, what was learned, what's next. See `CLAUDE.md` for the convention
and how this relates to `.claude/decisions/`.

---

## 2026-08-19 (evening) — manuscript-experiment redo pipeline

- Fixed two blockers found by the 100-cell manuscript runs: (1) phototroph
  respiration (both variants) drained bottom-bed NH4 negative — added
  depletion-protection Monods (NH4 1e-6, HPO4 1e-10) in both ports, suites
  green; (2) the negativity guard's diagnostic print crashed on index 0 —
  now reports cell/z/state-column/value.
- manuscriptExperiments driver: per-family detachment (biofilm 0.14 vs
  pathogen 1.4e-5·sqrt(v/18), the legacy srun split), Kappa=1e-7 (Table 2),
  cosmos account lu2026-2-100 (SSF; lu2025-7-124 is HDG — Jaime).
- Cosmos: first submission failed (module not found — sbatch needs
  `#!/bin/bash -l`); resubmitted as 3524225 (E3 summer 90 d, 500 cells) +
  3524226 (array E1/E2/E4/E5/E6/E10, afterok). Local 100-cell E3 canary
  running in parallel.
- sensitivity.tex updated to the v2 campaign numbers (corrected baseline,
  25 params, beta excluded).


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

## 2026-07-30 — sensitivity infrastructure repaired; Sobol restored to k=10; light sweep answered

Session ran overnight from 2026-07-29. Ten commits on `julia-port`, `e9e9a88..c3682e0`.

### What we accomplished

All items **ran and passed** unless marked otherwise.

- **`e9e9a88` Sobol flush + partial writes.** Job 3427791 ran 14.8 h with a 66-byte log and
  no output: Julia buffers stdout under Slurm redirection, and `sobol_indices.csv` was
  written only after the last parameter. Flush after each milestone with elapsed minutes;
  persist each QoI vector (`y_A`, `y_B`, `y_AB<i>`) as it completes. Partials go to
  `sobol_indices_partial.csv`, NOT `sobol_indices.csv`, which holds the complete Jul 22
  N=64 result — decision `.claude/decisions/2026-07-29-sobol-partial-file.md`. Verified with
  an N=2 run; Jul 22 file restored byte-identical.
- **`597fd12` Sobol resume.** Reuse persisted vectors when `resume` is passed; mature state
  built lazily so a full-cache resume costs nothing. Guarded: any change to
  N/k/params/ncells/tpost/nframes/tmature/seed aborts, because sampling changes make cached
  evaluations invalid — decision `2026-07-29-sobol-resume-guard.md`. Four checks passed:
  resume after deleting `y_AB4..8` gives byte-identical `sobol_indices.csv`; full-cache
  resume performs zero evaluations and zero mature builds; changed `ncells` aborts exit 1
  leaving results untouched.
- **`64636a8` light-profile step-5 aggregator** (`light_profile_report.jl`) plus the model
  and sweep drivers, until then untracked while running on cosmos. Refuses to report a
  direction from runs with `flag != OK`. Plan step 2 measured at last: **~1090 s per
  simulated day** at ncells=500, dt=3e-7 (2345 s cold), so cost was never the risk.
- **`64a2d9a` Sobol restored to k=10.** `velocity` and `beta_porosity` were excluded on a
  clogging rationale the cohesion fix made obsolete — all 27 parameters now report
  `flag = OK` in all three scenarios — while post-fix OAT ranks `beta_porosity` **first** in
  `pulse` (0.566) and `velocity` third (0.515). Jaime approved restoring both; `zeta_0`
  stays out at 0.045. Corner test: **0 of 20** range-extreme evaluations non-OK. β is
  sampled directly, not as the gap `1−β`, because log-uniform sampling of the gap is a
  different distribution over β and would break comparability with Morris. Decision
  `2026-07-29-sobol-restore-clog-params.md`.
- **`027f44f` o2pen censored + new o2dep.** Morris gave µ* exactly 0 for every parameter but
  `velocity`. Measured at ncells=30/pulse/1 day: O2 falls only to 60% of its maximum and
  0 of 62 cells reach the 1% threshold, so `findlast` returned the domain bottom as a
  constant. Now NaN when censored, plus `o2dep = 1 − min/max` as the uncensored companion
  (0.425/0.396/0.500 across the velocity range). Decision `2026-07-29-o2pen-censored.md`.
  **See the correction below — this diagnosis was incomplete.**
- **`2962a0c` slurm scripts tracked.** They existed only on cosmos, tracked nowhere, and two
  copies of `sobol-resume.sbatch` had already diverged unnoticed. Superseded 48-core draft
  removed; committed copy verified identical to the deployed one. README records the
  deployment direction and the two scheduling lessons.
- **`d91fcb6` cohesion discrimination tests.** Suite **216 → 243 passing**, all four golden
  suites still green. Asserts every preset's potential equals the published dψ/du at its own
  ζ₁ and differs from the double well at φ=0.01, and pins the stale-closure trap. The defect
  was wiring, not algebra.
- **`63047f8` light sweep result.** See below.
- **`b64273c` Morris flush + partial writes**, same two defects, same separate-file fix.
  Verified with a 2-trajectory run; Jaime's `mu_sigma.csv` restored md5-identical.
- **`c3682e0` Morris re-screen job.** 24 trajectories (not 6), 552 runs, `-c 2` because the
  driver is single-threaded.

### Light–biofilm profile study: answered

Job 3428150, all five amplitudes COMPLETED, `flag = OK`, 4:47–5:40 wall each.

**Jaime's hypothesis is contradicted.** Raising irradiance does not put more biofilm above
z = 0: supernatant mass falls 31% from its optimum, total mass 26%. There is a genuine
interior optimum at **A ≈ 0.05**, close to the A ≈ 0.018 predicted for `I_eff` = 1 at the
supernatant top, and both QoIs are non-monotonic across it — the photoinhibition signature,
measured.

**But QoI (a), the peak location, is the wrong observable.** `z_peak` is +2.00 mm at every
amplitude because φ_b nearly doubles across the sand surface (0.0151 → 0.0267 over one
cell) and then decays monotonically. There is no interior maximum; the argmax is the first
in-bed cell, pinned by the porosity ramp. So "the peak moves deeper" is **untested, not
refuted**. Use M(z<0) or a mass-weighted centroid instead. Numbers in
`.claude/plans/2026-07-29-light-biofilm-profile.md`.

### Corrections to claims made during this session

- **`o2pen` has two failure modes, and the committed fix addresses only one.** The diagnosis
  above (censoring: threshold never crossed, `findlast` pins to the last cell) was verified
  at ncells=30/pulse/tsim=1.0. But in the *Morris* regime the threshold **is** crossed
  inside the domain, and µ* is 0 because `o2pen` returns a **grid-quantised** depth `z[idx]`:
  unless a perturbation moves the crossing a whole cell (dz = 32.8 mm at ncells=30), every
  elementary effect is exactly 0. The running re-screen confirms it — `o2pen` comes back
  **0 NaN, 0 nonzero**, not the all-NaN I predicted. Quantisation is very likely the reason
  `velocity` was the only non-zero parameter in the Jul 21 file too. `o2dep` is unaffected
  and returns 21 of 22 non-zero. Quantisation is still unfixed.
- **The "rank agreement" claim in `reports/sensitivity-status-2026-07-29.typ` used the
  pre-fix Morris file and will need revising.** At 18/24 trajectories the post-fix Morris
  `Lmean` order is velocity 10.27 > **sand_pathogen 2.08 > dispersivity 1.85** >
  attach_sand 1.79, whereas pre-fix it was velocity > dispersivity > sand_pathogen >
  attach_sand. Ranks 2 and 3 swap, so post-fix Morris and post-fix OAT `startup` no longer
  agree on the full top four. Provisional until the run finishes.

### Cross-method validation emerging

Morris's σ/µ* flags are being confirmed by Sobol as the indices land:

| parameter | Morris σ/µ* (pre-fix) | Sobol S_i | Sobol S_Ti | reading |
|---|---|---|---|---|
| `attach_sand` | 2.22, strongly interacting | −0.012, CI straddles 0 | **0.296**, CI clear of 0 | almost entirely interaction — confirmed |
| `sand_pathogen` | 1.47, interacting | **0.379**, CI clear of 0 | 0.566 | real first-order plus interaction — confirmed |
| `dispersivity` | 0.76, near-additive | −0.001 | 0.003 | additive and negligible — consistent |

Also: at N=64 no first-order index resolved. At N=256 `sand_pathogen`'s S_i CI is
[0.218, 0.534], clear of zero. The larger N is buying what it was meant to buy.

### Plan for next session

1. **Read the two running jobs first** (see in-flight below). `sobolr` needs ~22 h more of
   its remaining 31 h; `morris` needs ~4.6 h of its remaining 10.3 h. Neither needs
   intervention — check `tail` of the logs and the `*_partial.csv` files.
2. **When Morris finishes, revise `reports/sensitivity-status-2026-07-29.typ`**: the rank
   agreement section, and add the post-fix magnitudes with 24 trajectories.
3. **Fix `o2pen` grid quantisation** — or retire the metric in favour of `o2dep`. A
   sub-grid interpolated crossing depth would resolve it; the current integer index cannot.
4. **Solver-level cohesion golden** — the outstanding piece. Measured: relative difference
   between the two potentials is 7e-17 from a seeded bump at φ_b ~ 5e-3, only 1e-10 with
   ζ₀ raised a millionfold, but **10–13% from a mature biofilm at φ_b ~ 0.16**. Cohesion
   flux scales with mobility u(1−u), so low φ_b is unconditionally blind. Needs a committed
   mature initial state anchored to MATLAB. **Direct seeding of
   `global_concentration.matrix` does not work** — the state is not self-consistent and
   trips the guards immediately; grow it with a pre-run.
5. **Re-examine the `--exclusive` trade** if Sobol needs to be faster in future. See risks.

### Open questions / risks

- **Sobol parallel efficiency is poor: ~1.2× on 16 threads.** 256 runs × ~52 s serial
  ≈ 222 min of work took 189 min wall. The node is shared (`--exclusive` was dropped to
  escape a 21 h queue wait, and Morris then landed on cn169 too), so memory-bandwidth
  contention is the likely cause. The job still fits its wall, so it was left alone, but
  **the exclusive/non-exclusive trade is the thing to revisit, not the thread count.**
- **`3427791` cost 14.8 h for nothing.** Cancelled along with `3428726`; both computed the
  superseded k=8 design. That loss is exactly what the partial-write fix prevents recurring.
- **The `slurm` echo line is stale**: `runs=$((N*10))` was hardcoded for k=8, so the k=10
  job prints 2560 where the truth is 3072. Log-only, cosmetic, unfixed.
- **ζ₁ mismatch still unresolved**: publication models use 0.005, presets 0.01,
  `results.tex` tabulates 1e-2.
- **`Logarithmic_OAT_Sensitivity_SSF_model.pdf`, `implicit-osmosis.md`,
  `worklog-2026-07-20.{md,html}`, `results/morris/`, and the two `MANIFEST.md` edits plus
  `zerocell_diag.jl` and `run_export_pathogen.m` are Jaime's uncommitted work** and were
  deliberately never staged.

## 2026-08-20 (late) — normalized light + Monod respiration; fleet reset

- Diagnosed I_opt units mismatch (absolute Wolf value vs dimensionless curves): whole
  supernatant photoinhibited, photosynthesis pinned to sand surface. Per Jaime:
  normalized_light=true (I(t) is already I/I_opt) and respiration light response
  changed to Monod K/(K+Î), K=1.0. See amendment 2026-08-20b in
  decisions/2026-08-18-phototroph-respiration-rate.md.
- Both ports committed/pushed; MPCSSF.jl synced; tests pass (light_inhibition asserts).
- Consistency reset: cancelled cosmos 3524363/64/3524433/34 (old regime), redeployed,
  resubmitted MATLAB 3524843/44 + Julia 3524845/46 (acct lu2026-2-100); local 100-cell
  suite + X2 wiped and relaunched. Two old-regime background runs (E10, X2) died from
  the wipe deleting their output dirs mid-run — expected, superseded.

## 2026-08-20 (midday) — kinetic-parameter audit

- Audited Table 2 growth rates + all 8 half-saturation constants against Wolf2007,
  Campos2006, Reichert2001 (now all in documents/): every value wrong vs its own
  citation; systematic "row slide" pattern. See decisions/2026-08-20-kinetic-parameter-audit.md.
- Probes: corrected rates alone change nothing (half-sat-limited); corrected half-sat
  set probing now (probeHalfSats.m); 2x2 season/respiration probes (probeSeason.m) running.
- ecomodel.tex temperature formula typo found: theta^(T/293-1) written, theta^(T-293)
  implemented (codes correct, manuscript wrong).

## 2026-08-20 (evening) — probe chain closed; co-author report written

- Probe chain: PG-excess seasons (winter 3.6x summer, structural), endogenous
  respiration 2x2 (deficit survives, dark growth gone), combined with corrected
  half-sats (net O2 producer, 13 mg/L diel swing), field influent 5e-4 from
  Campos2006b Fig 1(a) p.888 (clean 1.5 mg/L deficit, 0.3 mg/L swing, schmutzdecke
  back at sand surface). Influent algae in Table B.1 is ~15x field strength — the
  boundary condition was the other half of the O2 anomaly.
- Report: reports/kinetics-audit-2026-08.{typ,pdf} — audit tables, probe chain,
  endogenous-respiration proposal, 5 recommendations. Probe data archived in
  analysis/probes/data/.

## 2026-08-21 (early) — cosmos baselines timed out

- msE3 3524843 + jlE3 3524845 (500 cells, PG-excess normalized-light baseline)
  hit the 24 h wall limit with NO output written (dirs created, no frames/cache):
  a 90-day 500-cell run does not fit in --time=1-00:00:00. Nothing worth
  retrieving; regime is obsolete anyway (endogenous-respiration decisions
  pending). When the model of record settles: resubmit with a longer limit
  (partition max permitting) or split the 90-day run into chained checkpointed
  stages via the existing mature-cache mechanism.

## 2026-08-21 (continued)

### What we accomplished
- **Adversarial plan executed** (fixed-influent winter/summer flip evaluation;
  plan ~/.claude/plans/unified-juggling-widget.md, 2 Opus-critic rounds):
  - G-probe (5e985e8, matlab-claude): tau_P ~ 0.05-0.07 d both seasons ->
    P-routing (A) dead; **P closure violated ~45%** — HET death releases
    0.0209 P / 0.0653 N per unit vs 0.0141 / 0.0248 consumed (ran, flag OK).
  - H-probe (c568e2b): nutrient-closed death rows cut the defect 10x but
    worsen the flip (W/S 2.03 -> 2.51). Adopt as bug fix only (ran, flag OK).
  - Cardinal (CTMI) TemperatureResponse implemented both ports (5e985e8 MATLAB,
    89558ec Julia): opt-in, mu(20 C) anchor; MATLAB testCardinal PASSED, Julia
    283/283 PASSED, goldens bit-identical.
  - theta_growth sweep on cosmos (deb3b80; first run under the new
    cosmos-only probe policy, slurm/theta_sweep.sbatch): W/S 1.81/1.76/1.72
    at theta 1.066/1.09/1.12; winter theta=2.0 asymptote (job 3528329):
    **W/S floor = 1.65 measured with growth annihilated** -> E pre-empted,
    no growth-side intervention can flip; Bernard & Remond blocker moot.
  - **Budget-residual method exposed** (ef4e5d2): the "growth" residual is
    mu-independent (~0.05) -> all residual-derived mechanism numbers
    quarantined (decision file 2026-08-20-fixed-influent-flip-evaluation.md,
    addendum 2026-08-21b). Direct stock W/S measurements unaffected.
- Cosmos 500-cell E3 baselines timed out at 24 h, no output (journaled earlier
  today); queue cleaned.
- Memory: new feedback memory probe-jobs-on-cosmos (Jaime: never run
  simulation probes locally).

### Plan for next session
1. Fetch cosmos array **3528356** results (running at save time):
   task 0 = zero-biology mass check (dM vs supply-export -> quantifies the
   budget-residual bias), task 1 = summer theta=2.0 asymptote (true summer
   growth contribution). Commit mck_*/theta_summer_th2 data; also commit the
   already-fetched analysis/probes/data/theta_winter_th2.mat (untracked in
   the worktree).
2. Close the fixed-influent evaluation in the decision file with the mass-check
   numbers; fold G/H/theta findings + death-row nutrient violation into
   reports/kinetics-audit-2026-08.typ and recompile.
3. Bring the co-author package to Jaime: model-of-record stack (corrected
   rates/half-sats, PHO-only endogenous respiration, O2-neutral +
   nutrient-closed death rows, normalized light, seasonal influent, theta
   fixes incl. PHO growth 1.066 / HET death 1.08 / k_rb n/a) -> on approval,
   promote to presets both ports, anchor, rerun manuscript suite (local 100c
   + cosmos 500c with longer wall limits or chained 30-d stages).

### Open questions / risks
- Budget-residual bias source unknown until 3528356 lands (suspect phase-vs-
  bulk concentration convention in probe mass integrals; solver influent flux
  itself verified correct: porosity * (q/porosity) * c_in).
- All per-capita growth/loss mechanism numbers from probeBudget-family are
  quarantined; W/S stock ratios and O2 outputs stand.
- Death-row nutrient violation (P/N minting) is a manuscript-level stoich fix
  awaiting co-author sign-off alongside the rest of the audit.
- Main-checkout uncommitted files (MANIFEST.md edits, zerocell_diag.jl,
  run_export_pathogen.m, stray PDFs/worklogs) are pre-existing/Jaime's —
  untouched by this session.

## 2026-08-22 — Julia consolidated into SSF.jl; season90 probe on cosmos

### What we accomplished
- **Winter/summer, reframed by two field papers** (both new in `documents/`):
  - Bae2023 (full-scale SSF, Korea): schmutzdecke live cell counts are
    **winter 1.44 vs summer 1.25 x1e9 cells/g** (Table 2) — the model's W/S > 1
    is not contradicted by the field. Summer *activity* is higher (DOC removal
    43 vs 28%, UV254 59 vs 38%, Fig 4) while particle removal is flat (99%,
    Fig 3). Influent seasonality is 0.50-0.79 across every channel — **above**
    the model's 0.49 flip threshold, so on that site the model would not flip
    either. Table 2 is confounded by run time (200 d vs 330 d since scraping).
  - Bellamy1985: the temperature test is **fixed-influent** — Filter 6 chilled
    to 5/2 C beside Filter 1 at 17 C on the same raw water (Table 1), i.e. the
    E1/E2 protocol. Coliform removal 97 -> 87% (5 C) and 99.6 -> 92%, SPC
    99.9 -> 90% (2 C), while Giardia removal stays >99.9% at both. All
    activity, no standing-stock measurement.
  - Consequence: the manuscript's `results.tex:165` claim is about *growth*
    while `fig:seasons-results` plots *standing volume fraction*. Recomputing
    the existing probe pairs in phi_b gives W/S = 1.05-1.18 at the peak and
    ~1.0 at the sand surface, against 1.9-3.4 in PHO mass — the flip is largely
    an algal-compartment statement, not a biofilm one.
- **Transport "mass creation" is a probe artifact, not a solver bug.** The
  scheme conserves `sum(eps_i * c_i * dz)`: fluxes are porosity-weighted and
  divided by `porosityCenters` (`simulate.m:498-513`), and the sum telescopes
  exactly. Every probe integrates `sum(c)*dz` with no eps weight, against
  boundary fluxes already in bulk units (`q*c_in`). The implied bed fractions
  (PHO 33%, HET 89%) match where each organism sits. Julia mirrors the same
  conservative scheme (`julia/src/simulate.jl:521-533`), as does the legacy
  solver. NOT yet confirmed by a re-run: the saved .mat stores `phoM` summed,
  so the eps-weighted integral cannot be recomputed retroactively.
- **cosmos array 3530219** (6 tasks, `probeSeason90.m`): 90 d at manuscript
  resolution, three variants x two seasons, isolating the duration/influent/
  kinetics confounds and reporting phi_b and mass metrics both raw and
  eps-weighted. Running at journal time.
- **Julia consolidated into `code/1d/SSF.jl`** — decision file
  `2026-08-22-julia-consolidation.md`. `MPC-SSF/julia` retired after merging
  `julia-split` (SSF.jl 283/283, goldens bit-identical); package renamed
  MPCSSF -> SSF; `matlab-claude` merged (clean, 359 files); slurm/ split, with
  the Julia scripts' billing account corrected lu2025-7-124 -> lu2026-2-100.

### Plan for next session
1. Fetch 3530219 and compare `manuscript` vs `corrected` vs `field` at 90 d.
   The decisive cell is task 0/1: if the *published* configuration also gives
   winter > summer, the inversion predates every correction and
   `fig:seasons-results` does not follow from the committed model.
2. Re-run `probeMassCheck` with the eps-weighted integral to close the
   transport question (`dM == supply - export` to solver tolerance).
3. Seasonal removal arm: no probe measures removal, and `runMarker` computes log
   removals at a hardcoded 19 C (`manuscriptExperiments.m:336`) from a
   summer-only mature cache. A winter arm would test Bellamy/Bae directly and
   needs a `mature30_winter`.

### Open questions / risks
- W/S ratios are NOT automatically safe under the eps correction: the factor
  `1 + 1.5f` is distribution-dependent, and winter/summer PHO profiles differ.
- The seasonal influent (winter = 20% of summer, Campos2006b) is now the only
  leg the flip stands on, and Bae2023's 0.50-0.79 is well above threshold.
- Mauclaire2004 (doi 10.2166/aqua.2004.0009) is the sole support for the
  manuscript's summer claim and is not yet in `documents/`; its title is about
  clogging (a rate), not standing biomass.
- Cosmos trees are rsync copies still carrying `julia/`; stale but harmless.

## 2026-08-22 (continued)

### What we accomplished
- **slurm/ split** (`a493c5f` here, `2c4099b` in SSF.jl). Nine of thirteen batch
  scripts were Julia jobs living in what is now a MATLAB-only repo; they `cd`'d
  into `/home/jrman/MPC-SSF/julia` (deleted) and invoked `analysis/*.jl` that
  exist only in SSF.jl (verified file-by-file). Eight moved; three defects fixed
  in the move:
  - `-A lu2025-7-124` -> `lu2026-2-100`. lu2025-7-124 is the HDG project's.
    Both are valid for jrman on lu48 (`sacctmgr show assoc user=jrman`), so SSF
    work was billing the wrong project **silently**, not failing.
  - absolute log/working paths -> submit-directory-relative, matching the MATLAB
    scripts' existing convention.
  - `pkgtest` job name `mpcssf-pkgtest` -> `ssf-pkgtest`.
  `sobol.sbatch` was deleted by Jaime mid-task and treated as deliberate (it was
  the item-4 suggestion); removed rather than moved, its 48-thread lesson kept in
  the SSF.jl README prose. READMEs split; MATLAB deploy target corrected to
  `MPC-SSF-manuscript`, which is where those jobs have actually gone.
- **Decision + journal records** (`285d54a`): `.claude/decisions/`
  `2026-08-22-julia-consolidation.md`, per this repo's CLAUDE.md. Jaime's
  previously-unstaged 2026-08-21 journal section was committed with it.
- **Status report**: `reports/project-status-2026-08-22.{typ,pdf}`, compiled.

### Verification status (explicit)
- **Ran and passed**: SSF.jl `Pkg.test` 283/283 after the merge, goldens
  bit-identical (`RESULT: MATCH`); MPC-SSF/julia 283/283 before deletion;
  SSF.jl 276/276 pre-merge; `checkcode` clean on `probeSeason90.m`; typst
  compile of the status report.
- **Written, NOT run**: all eight moved slurm scripts. The account and path
  corrections have never been exercised by an actual submission — the first
  Julia job submitted from SSF.jl is the real test.
- **Written, NOT verified**: the eps-weighted mass integrals in
  `probeSeason90.m`. They are the proposed fix for the transport artifact but
  have not yet been checked against a case with a known answer.
- **In flight**: cosmos 3530219, 1 h 56 m into an 8 h wall, no output written.

### Plan for next session
1. Fetch cosmos **3530219** (6 tasks, `s90_*.mat`) into the matlab-claude
   worktree's `analysis/probes/data/`. Compare `manuscript` vs `corrected` vs
   `field` at 90 d in BOTH phi_b and PHO mass. Decisive cell is task 0/1: if the
   published configuration also gives winter > summer, the inversion predates
   every correction and `fig:seasons-results` does not follow from the committed
   model.
2. Re-run `probeMassCheck` with the eps-weighted integral; confirm
   `dM == supply - export` to solver tolerance and lift or confirm the
   budget-residual quarantine.
3. Recompute the W/S floor (1.65-1.88) eps-weighted — it is NOT automatically
   safe, the factor `1 + 1.5f` is bed-fraction dependent and the seasonal PHO
   profiles differ.
4. Seasonal removal arm (needs a `mature30_winter`): no probe measures removal,
   and `runMarker` is hardcoded to 19 C (`manuscriptExperiments.m:336`).

### Open questions / risks
- `probeSeason90.m` and `season90.sbatch` remain **untracked** in the
  matlab-claude worktree — the probe currently running on cosmos is not in git.
- Cosmos trees are rsync copies (no git) still carrying `julia/` and the old
  slurm scripts. Harmless now; redeploy before the next Julia run and let cosmos
  resolve its own Manifest (Julia <= 1.10.4 there vs 1.12.6 locally).
- `stash@{0}` in SSF.jl holds the superseded duplicate rename; safe to drop.
- Mauclaire2004 (doi 10.2166/aqua.2004.0009) still not in `documents/`; it is
  the sole support for the manuscript's summer claim and concerns clogging (a
  rate), not standing biomass.
- Nothing is pushed in either repo: SSF.jl is 5 commits ahead of origin/main,
  MPC-SSF julia-port is 48 ahead of origin/julia-port.
- Death-row nutrient violation (~45% P closure) still awaits co-author sign-off.

## 2026-08-22 (later) — Liebig-limitation diagnostic; cover-sweep plan

**Done.**
- Read all 13 PDFs in `SSF/documents/` and wrote `reports/literature-synthesis-2026-08.typ`
  (cross-paper discrepancies + model comparison + qualitative response matrix). Critic pass
  found four arithmetic errors of mine (three percent→log conversions, one bed-depth
  conversion) plus a wet-volume/dry-mass conflation in the φ_b comparison; all corrected
  inline with the challenges boxed in the report.
- Planned the referee-response experiments:
  `.claude/plans/2026-08-22-cover-sweep-and-referee-experiments.md`.
- Implemented `RecordLimitation` in `src/@State/simulate.m` (opt-in, default off).
  `evaluateReactions` now returns `[rx, monod, lim]`. Verified **bit-identical** three ways
  (pre-edit via `git stash`, post-edit flag off, post-edit flag on) across all 31
  concentration arrays.
- New probe `analysis/probes/probeLimitDiag.m`.

**Learned — two findings that reorder the plan.**
1. *The limiting nutrient was mis-assumed twice.* Manuscript stack: **NH₄**. Corrected stack:
   **IC** (supernatant biofilm) and **HPO₄** (flowing). Not the phosphorus-recycling story in
   `2026-08-20-winter-summer-flip.md`. The planned `+N` add-back arm would have relieved
   nothing on the corrected stack.
2. *Bed photosynthesis is absent, not small.* Î = 0.0608 at z = 0 → **1.39e-27** one cell in,
   on both stacks. `LightAttenuationCoeffSand = 1500` gives a ~5 mm euphotic depth against
   dz ≈ 9.95 mm. Campos2002's 0–2 cm sampling layer is two optically dark cells.

   Combined with min-Monod **0.703** in the top 0–2 cm on the corrected stack, the covered/
   uncovered null has **two regionally separate causes**: the supernatant is lit but
   nutrient-capped (≈0.005), the bed top is nutrient-replete but unlit. Neither factor alone
   explains it.

**Next.** Respecify Stage B (`+N` → `+IC`, add `+ALL`); promote Stage B′ (η_sand ÷ 10,
dz ÷ 3) to scheduled and run alongside B. Build `coverCases.m` / `probeCover.m` /
`collectCover.m`, then ship to cosmos (`-A lu2026-2-100`, `lu48`) — re-rsync first, the
remote trees are stale. Nothing beyond the smoke/diagnostic runs has been run locally.

### Cover-sweep infrastructure (same day, autonomous block)

Built `analysis/probes/{coverCases,coverTag,coverConfigFields,probeCover}.m`,
`analysis/collectCover.m`, `slurm/cover_sweep.sbatch`. Registry-driven: 29 unique cells
(A 5, B 16, B' 4, C 6), tagged by physics so the same configuration requested by two stages
runs once. Ladders and pairs share one reduction path. Full build log:
`SSF/reports/cover-sweep-build-2026-08-22.typ`.

**Two bugs the smoke test caught, neither visible to checkcode:**
1. `saferatio = @(a,b) a/b*(b>tol) + NaN*(b<=tol)` returns NaN for *every* input — `NaN*0`
   is `NaN`. Replaced with a branching `safediv`.
2. **The dark floor is on the wrong stack.** `modelLund.m` zeroes
   `MinimumLightFactor` only inside its `PhototrophRespiration > 0` branches (:119, :146).
   The manuscript stack passes that option and gets 0; the **corrected** stack appends
   `phoEndog` by hand, never enters the branch, and **keeps the 0.01 floor**. So a fully
   covered filter still photosynthesises at 1% of optimum — a fifth, unintended blocker on
   the very contrast the sweep measures, and backwards from what
   `results-redo-poc-2026-08.typ` records. `probeCover` now zeroes it explicitly.
   **`probeSeason90.m` has the same defect**, including the runs behind cosmos job 3530219.

**Verified:** lint clean ×5; 29/29 unique tags; solver bit-identical three ways; end-to-end
chain to 3 CSVs; stale-artefact assertion fires on a tampered cfg; `bash -n` and array
bounds cross-checked.

**Not run:** anything at scale. Needs an rsync to cosmos (remote trees stale) and Jaime's
ControlMaster socket. Open co-author question: the dark-floor defect affects season90 data
already on the cluster.

### Cosmos: recovered 3530219, deployed cover sweep, launched Stage A

**Recovered the "decisive cell".** `sacct` showed all six tasks of season90 array 3530219
FAILED (1:0) after ~3 h — but every task had printed `flag: "OK"`, `tFinal: 90` and saved
its `.mat`. Cause: **the same `fprintf(["a" "b"], ...)` string-array bug** as in probeCover,
at `probeSeason90.m:174`. The runs are sound; only the final summary line died. Fixed,
linted, redeployed. Two probes have now hit this independently — worth a lint rule.

All six fetched to the worktree's `analysis/probes/data/`.

**The seasonal answer (W/S at 90 d, 100 cells, fixed influent):**

| variant | phib_max | phib_int(eps) | PHO(eps) |
|---|---|---|---|
| manuscript | **1.125** | 1.673 | 3.432 |
| corrected  | **1.091** | 1.545 | 2.091 |
| field      | **1.027** | 1.266 | 1.955 |

Winter > summer in **every** variant including the manuscript arm. Per the criterion in
`project-status-2026-08-22.typ` §4, the inversion therefore predates every correction and
the published figure does not follow from the committed model. Caveats: the magnitude in the
plotted quantity is small (12.5% peak, 2.7% at field influent — inside the predicted 5-18%
and consistent with Bae's 1.15 +- noise), and the "manuscript" arm carries the 2026-08-20
normalised-light convention, so this does NOT isolate the light-convention change.

**New alarm only a long run could expose:** the corrected stack goes **near-anoxic** at 90 d
(bed min 0.024 / 0.041 mg/L) against 3.34 reported at 10 d. Elemo2024 measured filtrate DO
never below 3.0 and no anoxia. Separately, field influent **loses the O2 deficit** entirely
(summer effluent 9.25 vs 9.10 influent = net production). The redo's -3.7 mg/L used
PG-excess respiration, not endogenous, so the two configurations must be reconciled before
either number goes to a referee. The O2 claim in `literature-synthesis-2026-08.typ` §5.2 is
provisional until then.

**Deployed and launched.** 237 files rsynced to `~/MPC-SSF-manuscript/` (no `--delete`;
s90 data verified surviving 6/6). Registry verified on the cluster. Submitted
**job 3530596**, Stage A, `--array=1-5`, all five running concurrently on cn135.
Stages B/B'/C deliberately held: the plan's gates are human.

**Correction to the anoxia note above, after checking the profiles.** The bed does not sit
near-anoxic; it **oscillates into anoxia nightly**. Day-mean bed minimum on
`summer_corrected` is a comfortable 7.0 mg/L, but **32.1% of (cell, frame) samples in the
bed over the final day are below 1 mg/L** (winter 37.5%), instantaneous minimum 0.024. The
manuscript arm is 0% below 1 mg/L (min 3.27-3.54) and **field influent is 0%** (min
5.86-6.53). Physically coherent: pore velocity 18 m/d traverses the 1 m bed in ~80 min, so
the diel surface signal reaches the whole depth. So it is a corrected-kinetics + bloom-
influent pathology specifically, invisible to both a 10-day run and a daily-mean diagnostic.

Also: `o2_bed_min` in `probeSeason90.m` is misnamed -- it minimises over **all** z including
the supernatant, so two of the six "bed" minima are supernatant cells.
