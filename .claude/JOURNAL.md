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

## 2026-08-23/24 — kappa restored, implicit dispersion, and a working elimination

Long session. Three durable results, one open physical worry, and a running theme
worth carrying forward.

### Results that stand

**1. The Cahn-Hilliard kappa term was inert and is now fixed** in all three trees.
`getCahnHilliardMatrices` assembled every operator with `Rows <= n0`, leaving block
(2,1) empty, so `mu = Psi'(u^n)` exactly and kappa never reached mu. Verified
state-independently (block nnz 0, `max|muCH - psi'(u)| = 0` for random u), and
corroborated by the log-OAT report ranking kappa 27/27 with EXACTLY zero
sensitivity. Fix is `+ n0` on `diffusion.Rows`. MMS order 1.999; goldens re-exported;
`Pkg.test` 285/285. **Every previously published figure, manuscript Fig. 7 included,
was produced with no cohesion at all.**

**2. Implicit dispersion** (`ImplicitDispersion`, level A -- no phi_f lag) gives
**7.11x fewer steps / 5.06x wall at N=500**, measured three times independently,
for 6.1e-05 accuracy cost. It only pays at fine mesh (1.08x at N=100). Binding moves
from flowing-P/dispersion to enclosed-L, and at N=500 advection overtakes reactions
(0.53 vs 0.46) -- a crossover predicted at N ~ 407 from independent N=200 data.

**3. The component elimination WORKS**, via `CohesionScheme="matched"`: Solver A's
u-row rewritten as the exact flux Solver B applies. One-step agreement 8.3e-17;
reconstruction of enclosed water 3.67e-09 (2.46e-08 rel), 29,000x better than shin.
The residual is a stale `dt` in `rhsEnclosedWaterA` (:604 vs :984), not a floor.

Also: `IsUpwinded` default false -> true (centred convection has no maximum
principle; production impact 7.66e-04), and `CohesionScheme="bailo"` implemented and
validated against Bailo's own five gates in a standalone prototype.

### The open worry, and it is the important one

**E3 (kappa=1e-6, N=500, 20 d) shows biofilm migrating UP into the supernatant.**
Bed biomass peaks at day 10 and then falls 40%; by day 20, 57% of the biofilm is
ABOVE the sand, reaching 34 cm with phi_b = 0.14 at 30 cm. That is not a
schmutzdecke.

Suspected mechanism: `Psi(u) = u^4/4 - 0.5*zeta_1*u^3` has ONE minimum at
u = 1.5*zeta_1 = 0.015 and no high-concentration well, so `mu = Psi'(u) - kappa*Lap`
drives every cell above 0.015 downhill and the cohesive flux acts as POSITIVE
diffusion (Psi'' = 3u(u - zeta_1) > 0). It disperses the mat rather than holding it.

NOT yet attributable to the kappa fix: the Psi' half of that flux was always active
even when kappa was inert, and kappa sharpens interfaces, so it should oppose the
spread. Job **3534207** (kappa = 1e-6 vs 0, N=200, 20 d) settles it -- kappa = 0
reproduces the pre-fix scheme exactly.

Separately, **winter biomass 2.288 vs summer 1.545 at 20 d**: the kappa fix has NOT
restored the published summer > winter ordering.

### The theme: the instrumentation was less reliable than the code

Repeatedly, findings I reported with confidence turned out to be artefacts of my own
diagnostics, each corrected only by further measurement:

- the free-running diagnostic froze the phi_b-proportional detachment sink, removing
  the mechanism that bounds the solution -- I defended that design twice before
  testing it;
- the Bailo operators omitted the porosity weights Shin's carry, wrong by 2.5x at
  the z = 0 roughness layer, which made Bailo look WORSE than the scheme it improves;
- a sign error in the matched block (1,2);
- Bailo exporting a v_b built from a different mobility than it solved with;
- relative metrics with near-zero denominators reporting Inf, 4.3e12, 1.3e-04 as if
  meaningful -- three separate times;
- "drift saturates" read off an oscillation whose peaks were rising.

Hypotheses refuted along the way: that the reaction source caused the bound
violation; that convex splitting explained Bailo's excess drift; that the source
treatment drove the accumulated drift. Each was a reasonable guess killed by a cheap
test, which is the process working -- but the rate suggests testing instruments
against something independent BEFORE trusting them. `bailoCH1D.m` (the standalone
prototype, passing Bailo's five gates) is the one reference not downstream of this
machinery.

### Next

1. **3534207** (kappa 0 vs 1e-6) -- does the kappa fix cause the upward migration?
2. **3533845** (20 d drift) -- shin half valid; bailo half predates the v_b fix.
3. Fix the stale `dt` in `rhsEnclosedWaterA`; expect the reconstruction at roundoff.
4. Re-run cover sweep and legacy match: they predate BOTH the kappa fix and the
   upwind default.
5. The summer/winter inversion is a manuscript-level question, not a solver one.

## 2026-08-24 — Three-scheme Solver A drift

Armed `plotSolverADrift` for all three cohesion schemes and ran them as a Slurm
array on cosmos (job 3534837), plus shin/bailo locally for a machine control.

- **Fixed a blocker first**: `simulate.m` free-running diagnostic gated on
  `if ~useBailo`, so `matched` was advanced by the SHIN operator (`lhsPar`). A
  matched arm would have silently measured shin. See
  `.claude/decisions/2026-08-24-matched-free-running-branch.md`.
- **Fixed a latent plotter bug**: `D = D{k}` inside `for k = 1:2` clobbered the
  cell array, so the peaks overlay had never rendered and the last tile's
  `cellfun(@max, D)` would have thrown.
- **Split runner from plotter**: new `analysis/probes/probeSolverADrift.m` runs
  one scheme; `plotSolverADrift` is now the plot half and generalises to N
  schemes. That split is what allows the array.
- **Result overturns the hypothesis**: bailo did NOT subsume matched. Drift:
  matched 0.00% << shin 9.59% < bailo 14.17%. Matched's zero is tautological
  (its u-row is Solver B's update), which is what component-elimination wants but
  costs the independent check.
- **New concern**: the authoritative phi_b differs by 1.6% between bailo and
  shin/matched — a modelling-level disagreement, unresolved.

Next: decide the reference scheme on accuracy (MMS / refinement), not on drift;
then the secant splitting; then the zeta_1 sweep {0.01, 0.05, 0.15, 0.27}.

## 2026-08-24 (later) — Campos2002 / Demir2017 vs the model

Compared the field data against Manriquez2026's figures and the current code.

- **Data target:** biomass IN the sand, top 1.5-2 cm, monotone logistic to ~97 d,
  monotone hydraulic-conductivity decline. Campos2002 top-2 cm holds ~53% of top-10 cm
  biomass; Demir2017 "mainly the uppermost 1.5 cm".
- **Current code:** 96.9% of biomass in the SUPERNATANT at 90 d; bed peaks at 10 d then
  loses 8.8x; top 2 cm holds only 10-16% of bed biomass. Kozeny-Carman proxy on the bed
  rises then FALLS -- the model cannot clog.
- **Published figures:** Fig 9(c) (total mass) is reproduced and monotone, which is why
  the redistribution failure was never caught -- the diagnostic that was plotted passes.
  Fig 9(a)/(b) is not reproduced even at 20 d (published bed ~0.55 / sup ~0.15; ours
  0.278 / 0.325, inverted).
- **Three distinct problems, three levels:** (1) a regression, current code only --
  bisect kappa activation / IsUpwinded / normalised light; (2) too diffuse within the
  bed, published model too -- a cohesion-closure issue, and independent empirical
  support for raising zeta_1; (3) covered/uncovered light insensitivity, present in the
  PUBLISHED figures -- see `.claude/decisions/2026-08-24-light-insensitivity-covered-uncovered.md`.

Next: bisect the regression before any calibration; nothing downstream is meaningful
while the biofilm leaves the sand.

## 2026-08-25 — Deep study (cosmos 3537763)

10 arms, N=500 (one N=1000), 20 d chained in 10 d legs. Report:
`reports/deep-study-2026-08-25.typ`. Figure:
`analysis/results/figures/deepstudy_2026-08-25.png`.

**Infrastructure built first, because chaining was the prerequisite:**
- `probeRestart.m` — differences a chained run against a continuous one, field by
  field. PASS: no snapshot field zeroed, worst difference 4.36e-05. Seeding the
  continuation dt (rather than cold-starting at 1e-8) improved it 3x from 1.26e-04,
  confirming the residue is the adaptive stepper re-ramping, not lost state.
- `stateFromFrame.m` / `rehomeState.m` — extracted from manuscriptExperiments'
  locals. `grab` silently substituted ZEROS for scalar-stored components; it now
  warns, and the snapshot carries an `allZero` list that probeChain treats as fatal.
- `probeChain.m` — runs in <=10 d legs, saves per leg, resumes by skipping completed
  legs. Saves profile history, bed/supernatant trajectories AND effluent for all 9
  species.
- `checkRunFlag.m` wired into every probe (see 2026-08-24).

**Results:**
- NUMERICAL: converged at N=500/MaxDt=3e-6. Bed biomass agrees to 4 sig figs across
  MaxDt 1e-6 -> 5e-5 (50x), and 0.74% between N=500 and N=1000.
- MECHANISM: the bed peak-and-decline is driven by the SUPERNATANT LAYER, not
  detachment. Decline vs final supernatant fraction correlates at r = 0.91. Doubling
  detachment gives the SMALLEST decline (1.09x vs 1.25x baseline) because it
  suppresses the layer. Detachment hypothesis (Leg A) refuted as primary.
- With DetachScale = 0 the run CLOGS at t = 6.46 d in cell 1, the top of the water
  column, at 99% biofilm. Detachment is the only bound on total biomass.
- VALIDATION: effluent O2 is invariant across the entire sweep (44.4-44.7% removal),
  matching Campos2002's covered/uncovered TOC/DOC result (25/23 vs 23/23% despite
  4x biomass). Corollary: effluent data cannot discriminate zeta_1, detachment or BC.

**Bugs found and fixed today:** leg-2 time offset in probeChain (simulate returns
ABSOLUTE frame times when the state carries a start time; adding tStart
double-counted, putting leg 2 at t=20..30).

Next: an arm with biofilm light-attenuation disabled, to split "shading" from
"substrate capture" within the supernatant-layer mechanism.

## 2026-08-25 (later) — ζ₀=5 to 40 d, N=1000 check, critique + adversarial review + code review

- **ζ₀ = 5, N=500 extended to 40 d** (legs 5–6, hourly frames). Does not clog: φ_b(0)
  saturates at 0.845. Bed keeps declining (1.12× from the 16.8 d peak); total biomass
  at steady state by 40 d; supernatant 41.5 %. Figure
  `analysis/results/figures/zeta0_5_n500_40d.png`. Record:
  `.claude/decisions/2026-08-25-zeta0-tradeoff.md` (the ζ₀ arms had no durable record
  until now — critic finding).
- **N=1000 at ζ₀=5, 10 d**: bed +0.37 % vs N=500, profile L2 0.25 %. Converged.
- **Mass closure CLOSED**: new `analysis/probes/probeMassClosure.m` (zero biology,
  ε-weighted stocks): residual 1.1e-9 of supply for HET and PHO. The 08-21c "transport
  creates mass" was the unweighted probe integral. Quarantine lifted for ε-weighted stocks.
- **Head-loss diagnostic** `analysis/headlossKozenyCarman.py` (post-processing only).
  H/H₀ crosses 3× at 4.5–5.0 d in EVERY ζ₀ arm, then falls as the bed loses biomass.
  Operational clogging is 10–20× too fast vs Demir2017 and non-monotone; ζ₀ irrelevant.
  Table: `analysis/results/headloss_zeta0_2026-08-25.txt`.
- **`.claude/CRITIQUE.md`** written (scoreboard vs goal, unverified foundations,
  ζ₀ honestly, light, OAT, process, proposals), then adversarially reviewed by a critic
  agent (§9) and answered (§11). Biggest hits taken: Bae2023 row was inverted (field
  W/S 1.15 SUPPORTS winter ≥ summer — the published summer ≫ winter figure may be the
  wrong target); OAT ran at 30 cells, not 100, and already includes ζ₀ and β; the
  "same root cause" link between covered/uncovered and winter/summer is a hypothesis,
  not a finding; `coverCases.m` already encodes the proposed cover probes (Stages B/B′/C).
- **`.claude/CODE-REVIEW-2026-08-25.md`** by a third agent (src/ only, report only).
  Should-fix: `Model.m:41` detachment arity; `pathogenModel.m:72` 1.4e-5 vs paper 0.14;
  dead paths `Model(Preset=)`, `State.VolumeFractions`, `Results.concatenate`; Solver A
  uses dt_old in adaptive mode; `SolverOptions.override` skips validation; CFL bound
  hard-codes the Lund layout; `modelLund.m:170` sqrt without abs. Hot-loop bloat list.
  Nothing applied.
- mg→µg Campos correction swept into `decisions/2026-08-24-light-insensitivity-*.md`.

Next (from CRITIQUE §7.4 as amended): run `slurm/bisect.sbatch`; add the
seasonal-influent arm to `probeSeason90.m` and run the pair; run `coverCases` Stages
A/B/B′/C at ζ₀=5; add a finite-value check to `checkRunFlag`; decide the code-review
should-fix list with Jaime; then the 90 d steady figure and the OAT at N=500.

## 2026-08-25 (evening) — photosynthesis was off; half-saturations corrected

- Field-influent 2×2 (`chain_z0w_5_n500_field{,_dark}`) and the dark manuscript-influent
  arm: lit = dark to 3 s.f. at BOTH influents. Root cause: `modelLund` K_NH4,PHO = 1.2e-2
  → Monod 1.7e-4 at the surface, growth 5.4e-4/d vs death 0.37/d. The 08-20 audit had
  found it; the preset never received the fix. Every run this week had photosynthesis off.
- **`modelLund.m` half-saturations replaced with the audited set** (HET O2 2e-4, DOM 4e-3,
  NH4 1e-6, HPO4 2e-5; PHO IC 1.2e-3, NH4 2e-5, HPO4 2e-5; hydrolysis 0.1). Old values
  kept in a comment. `testRespiration`/`testCardinal` pass; smoke simulate OK. Decision:
  `.claude/decisions/2026-08-25-half-saturations-corrected.md`. NOT changed: MarkerGrowth
  (`pathogenModel.m:55`, drives the OAT), respiration O2 half-sat, SSF.jl mirror + goldens.
- `probeChain` gained `LightScale`, `Influent`, and saves `*_ABORTED.mat` on non-OK flags.
- ζ₀=1 N=500 re-run: CLOGGED 17.00 d as a one-cell spike at z=0 (0.986; first bed cell
  0.648). N=1000 arm running.
- Influent check: PHO 1e-2 is ~15–20× field (Campos2006b), HET 2.68e-3 rests on a 6 µm³
  cell; IC derivation inverted in the text. Field set: HET 1.5e-4, PHO 5e-4, HPO4 5e-6.

ALL pre-2026-08-25-evening results are on the crippled phototroph and must be re-run.

## 2026-08-25 (night) — Reichert respiration, kinetics audit, preset fixed

- Corrected-K 2×2: photosynthesis ran (5 g/m³ O₂ diel, noon 10.1 > influent) but biomass
  unchanged; enclosed NH₄ 70× below flowing. Diagnosed asymmetric O₂ stoichiometry
  (growth +0.93, death +0.20, PG-excess sink a biomass source).
- **`RespirationForm="reichert"`** added to `modelLund` (RWQM1 (10): PHO −1, O₂ −0.9301,
  IC +0.36, NH₄ +0.06, HPO₄ +0.01; 0.1/d, θ 1.047, K_O₂ 2e-4). Decision
  `2026-08-25-reichert-respiration.md`. `probeChain` uses it.
- Six Reichert arms (`z0w_5_rr_*`): phantom carbon gone (bed biomass halves, POM export
  13.8→5.7 g/m³ < particulate import); **first light effect on biomass** — manuscript-influent
  lit mat 2.6× dark; but O₂ still supersaturated (mean 11.9) because heterotrophs were off
  (μ_HET 0.018/d uncorrected), and at field influent lit ≡ dark (+N changes nothing: 91 %
  of phototrophs sit in the dark sand). Figure `zeta0_5_reichert_influent_x_light.png`.
- ζ₀=1 Dirichlet vs Neumann (corrected K): clog 15.18 vs 14.07 d; BC acts on one cell via
  κ/dz² = 0.25; the 0.84→0.66 jump at z=0 is Solver-B matrix transport (wall), forms at
  day 5–6 when the supernatant catches the saturated bed (0.69). Figures
  `z0_1_neumann_vs_dirichlet{,_t468}.png`.
- **Kinetics audit vs Wolf2007 Table VI + Reichert2001 Table 8, preset fixed**
  (`2026-08-25-kinetics-wolf-reichert.md`): μ_HET 0.018→2.0 (θ 1.0725), μ_PHO 5.5→2.0,
  d_HET 2.0→0.4 (θ 1.0725), d_PHO 0.4→0.1 (θ 1.047), k_hyd 0.09→3.0 (θ 1.0725); half-sats
  re-checked (Wolf/Reichert disagree on N: HET 10⁵×, PHO 6×). Tests pass, smoke OK
  (0.5 d: O₂ max 9.12 ≈ influent, HET matrix 0.2 already).
- Influent check: Chan2018 PDF holds no NH₄/NO₃ numbers (Table S1 not in documents/);
  nitrate absent from the model. Plan `2026-08-25-covered-uncovered-contrast.md`.

ALL runs to date predate the rate fix; the 2×2 must be re-run on the final preset.

## 2026-08-26 — Campos route: respiration merged into death

- Final-preset 2×2 + respiration sweep (`z0w_5_fp_*`): O₂ flipped from supersaturated
  to near-anoxic at manuscript influent (0.3–1 g/m³; bed mineralises 12.7 g/m³ of
  influent particles), 8.5 at field influent; POM export 14 → 0.09 g/m³; HET now 28 %
  of the biofilm; lit mat 2.9× dark at manuscript influent. Sweep r = 0/0.1/0.276/0.5/1.0
  flat and non-monotonic: respiration is not the O₂ lever. Figure
  `zeta0_5_finalpreset_2x2_respsweep.png`.
- DOM 1 mg/L at field influent: 2.7 % removal at 20 d — K_DOM 4 mg/L freezes HET growth
  at its death rate (0.37 vs 0.37 /d). Lit ≡ dark.
- Campos2006/2006b read: covered/uncovered contrast is imposed structurally (schmutzdecke
  layer + 30 % carbon export to top 2 cm), uncovered bed biomass peaks then declines as the
  mat takes over filtration, thicker mat → lower head loss (Nakamoto).
- **Preset: d_PHO = 0.276 /d, θ 1.08 (Campos k_ra), no respiration reaction**
  (`2026-08-26-campos-route-no-respiration.md`). Tests pass, smoke OK.

## 2026-08-26 (early) — sweep wave 1: the lock is the flowing→enclosed transfer

- Partner agent (Opus 5) proposal `.claude/PARTNER-2026-08-26.md`: target table T1–T9 with
  sources, 22-arm sweep, verdict on walking back O₂ modifications (no; the O₂ problem
  inverted), risks. Its live trap fixed: `modelLund` MinimumLightFactor default 0.01 → 0.0
  (`testRespiration` updated).
- N=200 vs N=500 (`*_n200` twins): bed/O₂/DOC/head loss within 2 % or 0.06 g/m³; the mat is
  30 % small and 2 d late. N=200 for bed/O₂/DOC screening; N=500 for anything about the mat.
- **K_DOM ladder** (field + DOM 1 mg/L): removal 2.7 → 8.5 → 10.2 → 10.8 % for 4e-3 → 1e-4;
  saturates. Enclosed DOM is 500× below flowing, enclosed O₂ 20× — transfer-limited.
- **Gate** (PHO_in 5e-5, DOM 3 mg/L, η_sand 150, N=500): nothing grows lit or dark
  (bed 0.00038 vs 0.00034, no mat). The model cannot grow a biofilm; all past mats were trapped.
- **Transfer ladder** (liquid TransportRate ×1/×10/×100/×1000 at field+DOM1, K_DOM 3e-4):
  DOC removal 10 / 74 / 95 / 95 %; effluent O₂ 8.29 / 7.28 / 6.93 / 6.93 (Elemo window 2–5
  consumed: ×10–×100 inside); bed +48–65 %; HET share 54 → 72 %. Saturates at ×100; ×1000
  costs 5× wall for nothing. Lit ≡ dark at every rung (no mat at field loading).
  600/d ≈ 400 µm diffusion length; ×10 ≈ 130 µm. With DOM_in = BDOC (~30 % of a 3–4 mg/L
  DOC), 74 % BDOC removal ≈ 22–25 % total-DOC removal = Campos2002.
- probeChain options added: KDOM, EtaSand, TransferScale (all recorded, guarded).
- Launched wave 2: PHO_in ∈ {5e-4, 1e-3, 3e-3, 1e-2} × lit/dark at transfer ×10, K_DOM 3e-4,
  DOM 1e-3, N=200.
- **Wave 2, load ladder** (transfer ×10, K_DOM 3e-4, DOM 1e-3, HET_in = 0.268·PHO_in, N=200,
  20 d; `chain_sw_P{5e4,1e3,3e3,1e2}_{lit,dark}`):

  | PHO_in | O₂ out (consumed) | H/H₀ 20 d | top-2 cm µg C/g dry / wet | mat | lit/dark |
  |---|---|---|---|---|---|
  | 5e-4 | 7.32 (1.8) | 1.04 | 105 / 26 | none | ≡ |
  | 1e-3 | 6.45 (2.65) | 1.08 | 187 / 47 | none | ≡ |
  | 3e-3 | 3.85 (5.25, min 3.68) | 1.28 | 398 / 99 | none | +2 % |
  | 1e-2 | 0.03 (anoxic) | 2.39, t(2×) 5 d | 701 / 175 | 0.0125 lit / 0.0056 dark | 2.2× mat |

  T1 (Elemo 2–5 consumed, never < 3) is met for PHO_in 1e-3–3e-3; T6 (H/H₀ 1.3–1.8 at 20 d)
  approached only at 3e-3; Campos's 60 µg C/g sits at ~5e-4 (dry basis) or ~1.5e-3 (wet).
  No supernatant mat forms below 1e-2 and lit ≡ dark there — the mat is a trapping threshold,
  not growth. Partner corrected my gate reading (it ran at transfer ×1); gate re-run at ×100
  (`g2_*`) and a phototroph-only arm (`g3_*`, HET_in = DOM_in = 0) running at N=500.
  F4–F6 (K_DOM 1.9e-3 / 4e-3 at ×10, 4e-3 at ×100; PHO_in 1e-3) launched.
- Partner §G: state variable is wet mass, RWQM1 stoichiometry is per dry OM — a ~4× yield
  error; fix = ρ_P/f_dry with influents scaled, φ_b path invariant. Needs Jaime's f_dry.
- **F4–F6** (PHO_in 1e-3, DOM 1e-3, 20 d): BDOC removal 86 % (K 3e-4, ×10) / 54 % (K 1.9e-3
  Campos, ×10) / 14 % (K 4e-3 Wolf, ×10) / 21 % (K 4e-3, ×100). As total-DOC removal (BDOC =
  30 % of DOC): 26 / 16 / 4 / 6 %. Wolf's ASM K is too high at any transfer; Campos's ksCd
  under-removes; Campos2002's 23–25 % sits at K_DOM ≈ 1e-3 with transfer ×10 (L_f ≈ 200 µm).
- **Gate at transfer ×100** (`g2_*`, N=500): identical to ×1 — nothing grows (bed 0.00038,
  no mat, enclosed O₂ 0.14 vs 8.96 flowing at z=0). Cause is the transfer FORM:
  transL ∝ φ_enclosed/β, so a thin film cannot be fed whatever the coefficient. This is the
  seeding lock. **Phototroph-only arms** (`g3_*`, HET_in = DOM_in = 0) aborted at 7 d with
  NaN in HET matrix at z=0: hydrolysis quotient POM/HET = 0/0 when HET ≡ 0 (a guard is
  missing in lookupQuotients). Re-launched as `g3b_*` with HET_in = 1e-7. Partner asked for
  an area-scaled transfer form and the operating point.
- **CORRECTION (08-26 morning).** `Frames.Concentrations` stores BULK values (per total
  volume) for every phase; the enclosed-phase *concentration* is bulk/φ_e. Every
  "enclosed/flowing" ratio I reported overnight compared bulk enclosed to bulk flowing, i.e.
  was low by the factor φ_e (0.01–0.08). Recomputed locally: kd_1e4 (×1) DOM ratio 0.03,
  O₂ 0.6 — still transfer-starved but 30×, not 500×; xf_100 ratio 0.94 (equilibrated);
  **g2 gate ratio 1.00 — the thin film was fully fed and still did not grow.** The
  "transfer form is the seeding lock" conclusion (and the partner's §H built on my numbers)
  is therefore WRONG as an explanation of the gate. The film form (D/L_f², now implemented
  behind `SolverOptions.TransferForm="film"`, exact relaxation, CFL term dropped) remains
  the physically justified replacement for the uncited 600/d, but it is not why nothing
  grows at low load.
- Local-concentration analysis of the gates: at the sand surface at noon, enclosed HPO₄ =
  5 µg/L → Monod 0.20 (binding), gross growth 0.33/d vs loss 0.256/d — net negative on the
  daily mean. Autotrophic growth from a clean start is P-limited at the 5 µg/L set from
  Chan2018 "<10 µg/L total P"; with trapped biomass recycled P reaches ~20 µg/L and noon
  net is +0.6/d (a bistability signature — Campos2006b's lag phase). Partner §K–N: no paper
  in documents/ reports an influent phosphate; our K_HPO₄ 2e-5 is Reichert's, while Campos's
  own ksp range is 1e-6–5e-5. Launched N1 (KHPO4 ∈ {2e-5, 5e-6, 1e-6} lit/dark) and N4
  (HPO₄_in 2e-5 scenario, lit/dark) at the operating point (PHO 1e-3, DOM 1e-3, K_DOM 3e-4,
  transfer ×10, N=200). Film-form A/B (`fa_film_lit`) running. `probeChain` gained KHPO4.
- **Film-form A/B** (`fa_film_lit` vs `sw_P1e3_lit`, PHO 1e-3, K_DOM 3e-4): film gives
  enclosed/flowing (local) 0.96–1.00 vs 0.13–0.25 at constant ×10; BDOC removal 95 vs 86 %;
  O₂ 6.29 vs 6.45; bed +3 %; and 30 % less wall (CFL transfer term gone). Film ≡ the ×100
  saturation, as predicted. With film, K_DOM must rise toward Campos's 1.9e-3 to land 23–25 %.
- Effluent particulates: HET 6-log, PHO 4-log removal vs Bae2023's 98.7–99.3 % (2-log) —
  attachment 547/d vs detachment 0.14/d is 10⁴ too retentive. `probeChain` gained
  `DetachForm` ("sqrt" | "linear", equal at 18 m/d) for a detachment sweep.
- **P half-sat ladder (N1) + P-replete scenario (N4)**: nothing. lit/dark bed ratio 1.004 →
  1.011 for K_P 2e-5 → 1e-6; HPO₄_in 2e-5 gives 1.009. φ_b(0) lit 0.200 vs dark 0.176 (+14 %)
  is the only light signal. Phosphorus is NOT the lock. Daily-mean arithmetic (Steele mean
  0.435 at the sand, 17 h of light): preset gross 0.42/d vs loss 0.256 → +0.16/d with P
  unlimited — growth is positive, yet no contrast: light acts on ONE cell (z=0; e-folding
  1.1 mm below, 3 h flushing above with attachment ∝ φ_b in the water column), so the
  phototrophs are retained only in the dark sand. The lock is retention in the lit zone.
  Launched growth/loss ladder on the film form: d_PHO 0.1, μ_PHO 3.0, both, both+η_sand 150,
  lit/dark (`gl_*`). `probeChain` gained MuPHO, DPHO.
- **Detachment sweep** (film, PHO 1e-3, 20 d): total particulate removal (HET+PHO+POM out vs
  in) is **98.5 %** at ×1 — already inside Bae2023's 98.7–99.3 %; the "6-log" was HET alone.
  ×3 → 95.2 %, ×10 → 63 %. Linear ≡ sqrt at this loading (φ_b max 0.19; v_f barely rises);
  the forms only diverge near pore closure. Top-2 cm carbon 217 → 77 (×3) → 29 (×10) µg C/g
  (dry basis); ×3 lands on Campos's 60 at the cost of 3 points of particle removal. Keep ×1;
  revisit linear at ζ₀ = 1 (interface spike) only.
- **Film transfer rejected as default** (`2026-08-26-film-transfer-rejected.md`): it moves
  the φ_b/HET maximum to 2–5 cm depth (constant ×100/×1000 do the same); Jaime: no
  reasonable simulation has φ_b peaking inside the bed. Cause: leaked hydrolysate advected
  downstream. Constant ×10 keeps the surface maximum. Film-form δ/μ arms killed and
  relaunched on constant ×10; two film diagnostics (N=500, DOM_in=0) left running.
- Growth/loss ladder (film): d_PHO 0.1 → bed +42 % symmetric lit/dark; μ 3.0 nothing;
  only η_sand 150 gives lit/dark 1.16 (bed), 1.08 (top 2 cm) and a first supernatant PHO
  in the lit arm. Retention in the lit zone is the lock.
- **Grid/light finding (decisive for the lit≡dark null).** Noon Steele factor per cell:
  N=200, δ=5 mm, η_sand=1500 — the sweep grid — has NO cell inside the roughness layer
  (the −0.5 cm cell has ε=1, no bare-sand attachment) and the z=0 cell is at L=0.16 (daily
  mean ~0.07), z=+0.5 cm at 0.00. There is no lit cell that can hold biomass; every N=200
  lit/dark null this morning is a grid artefact of the light–roughness coupling. N=500 has two
  roughness cells (−0.4, −0.2 cm; ε 0.88/0.64; L 0.85/0.54) → the N=500 field arms showed
  lit/dark 1.02. η_sand=150 lights z=0..+1 cm (L 0.79→0.42) → lit/dark 1.16. δ=2 cm at
  N=200 gives a lit ramp cell (−1.5 cm, L 0.64, ε 0.85) but lit ≡ dark (1.000) because
  growth there is P-limited (influent P, no recycling). Nothing has yet combined a lit
  attaching zone with P relief. Launched `lz_{n500,d2cm,eta150}_kp_{lit,dark}` (K_P 1e-6).
- δ=2 cm (constant ×10): φ_b(0) 0.385 (2× δ=5 mm), first supernatant PHO 0.0014, lit≡dark.
  mu ×1 pair: lit/dark 1.004. mu ×10/×20 arms slow (stiff), still on leg 1 after 30 min.
- mu ×10/×20 arms killed at 39 min (leg 1 unfinished, stiff); not part of the deliverable.
- **Lit-zone × P arms (K_P 1e-6, constant ×10, PHO 1e-3, N=200)**: δ=2 cm → lit/dark 1.001
  (nothing, even with a lit attaching ramp cell and P relieved); η_sand=150 → 1.133 bed,
  1.131 top-2 cm. P relief adds nothing over η_sand alone (1.16 on film). Light into the
  sand is the only lever and gives 13–16 %, not Campos's 4×. Structural: no mechanism holds
  phototrophs in the light. N=500 pair pending. Deliverable: `reports/results-2026-08-26.pdf`.
- **N=500 lit-zone × P pair**: lit/dark bed 1.036, top-2 cm 1.051, phototrophs above the sand
  2.4× (0.00245 vs 0.00101 kg/m²), diel O₂ max 6.57 vs 6.48 — the first light response at
  realistic loading, carried by the two resolved roughness cells; the bed itself moves 4–5 %.
  Report updated and compiled: `reports/results-2026-08-26.pdf`. No jobs running.
- **Parameter set chosen** (two independent picks, `.claude/PICK-A/B-2026-08-26.md`, reconciled
  in `reports/results-2026-08-26.pdf` §8): PHO_in 3e-3 / HET 8e-4 / DOM 1e-3 / HPO₄ 5e-6;
  K_DOM 3e-4 (1e-3 recommended), K_P 1e-6; μ_PHO 2.0, d_PHO 0.276, no respiration; transfer
  constant ×10; detachment ×1 sqrt; η_sand 1500, δ 5 mm; N=500; ζ₀ 5, ζ₁ 0.27. Behaviours:
  O₂ met (5.3 consumed, min 3.68), lit/dark met-small (N=500: 1.036/1.051/2.4× above sand),
  monotone met, steady state partly, H/H₀ 1.28 partly, surface maximum met, DOC met/partly,
  particulates met. Next: this set lit/dark at N=500 to 104 d on cosmos.
- **ζ₀ = 1 with LINEAR detachment** (chosen set, N=500, lit/dark; `z01_set_lin_*`): no clog.
  φ_b(0) 0.514 (interface cell, neighbours 0.38/0.32), supernatant 0.3 % of biomass (vs
  30–40 % at ζ₀=5, manuscript influent), bed monotone to 20 d, H/H₀ 1.29 (in band),
  O₂ 3.88 mean / 3.58 min, top-2 cm 352 µg C/g dry, lit/dark bed 1.045, top-2 cm 1.053,
  PHO above sand 2.46×. Linear detachment caps the z=0 pile-up (k_det 0.9–2.8/d there vs
  0.35–0.65 with sqrt). This supersedes ζ₀=5/sqrt as the working set.
- CRITIQUE §12 added: buoyant mat export at the top boundary is a missing mechanism.
  Detachment moves biomass into the flowing phase and it re-attaches ~3 cm into the bed; the
  model has no upward route and no overflow. Nakamoto's mat is lifted by photosynthetic O₂
  bubbles and drains out through the overflow with the solids it trapped — which is why his
  uncovered filters clog LESS while Campos2002's clog EARLIER (interstitial biomass). We are
  structurally locked into the Campos2002 regime; the manuscript must claim only that.
- **E4: 2× field influent to 60 d** (`fld2x_{lit,dark}`, 6 legs, 1.71 h per arm). First run to
  reach steady state: last 10 d change +0.14 % bed, −0.006 g/m³ O₂. O₂ consumption 2.96 mg/L
  (min 6.07) — Elemo window, at steady state; BDOC 89 % ≈ 27 % total DOC — Campos band; bed and
  head loss monotone; surface maximum; supernatant 0 %. **H/H₀ plateaus at 1.10 vs Demir's
  2.1–4.2 at 55 d — measured, not extrapolated: this set never clogs.** Light contrast keeps
  growing while all else saturates (PHO above sand 2.54×, bed 1.040) — evidence the light
  response is growth, not trapping. Recorded as E4 in EXPERIMENTS.md; figure
  `analysis/results/figures/fld2x_60d.png`.
- **E5: μ_PHO ×10/×20 at N=500** (`mu500_x{10,20}_{lit,dark}`, ζ₀=5, √-detachment — launched
  pre-E3). A 20× growth rate gives lit/dark bed 1.083 vs 1.025 at the preset, and 40/d is
  identical to 20/d — saturated. Supernatant stays 0.5 %; H/H₀ and O₂ unmoved. μ_PHO is closed
  as a lever: the contrast is limited by the ~4 mm of lit attaching zone, not by growth rate.
  Cost: 2.8 h and 5.0 h per lit arm (5–9× the preset) — stiff reaction term. No figures.
  Recorded as E5.

## 2026-08-27

### What we accomplished

Branch `julia-port`; **no commits made this session** — all work is uncommitted (13 modified
tracked files, 112 untracked). Last commit remains `88c8b15`.

**Model corrections (all with decision files, all uncommitted):**
- `src/presets/modelLund.m` — half-saturations replaced with the Wolf2007/Reichert2001 set
  (`2026-08-25-half-saturations-corrected.md`); growth/loss rates and temperature factors
  audited and replaced (`2026-08-25-kinetics-wolf-reichert.md`: µ_HET 0.018→2.0, µ_PHO 5.5→2.0,
  d_HET 2.0→0.4, k_hyd 0.09→3.0, θ from Reichert β); phototroph respiration merged into a single
  Campos k_ra = 0.276/d loss with no respiration reaction
  (`2026-08-26-campos-route-no-respiration.md`); dark-growth floor default 0.01→0.
- `src/@State/simulate.m` — film-form liquid transfer added behind
  `SolverOptions.TransferForm="film"` and then **rejected as default**
  (`2026-08-26-film-transfer-rejected.md`: it moves the φ_b maximum 2 cm into the bed).
- `analysis/probes/probeChain.m` — options added: `LightScale`, `Influent`, `Respiration`,
  `KDOM`, `KHPO4`, `EtaSand`, `TransferScale`, `TransferForm`, `DetachForm`, `MuPHO`, `MuHET`,
  `DPHO`, `Delta`; all recorded in `rec` and covered by the resume guard; partial trajectories
  now saved as `*_ABORTED.mat` on a non-OK flag.
- `analysis/headlossKozenyCarman.py` (new) — Kozeny–Carman bed head-loss diagnostic.
- `analysis/testRespiration.m` — updated for the zero dark floor. Tests: `testRespiration`
  and `testCardinal` **ran and passed** after every preset change; `checkcode` clean.

**Experiments** (all in `EXPERIMENTS.md`, E1–E6, with commands, parameters, wall times,
output files and figures):
- E1 field / 2× field at the working set, 10 d — O₂, DOC, particulates in band.
- E2 µ_HET ladder — **rejected as a lever**: lowering µ_HET buys O₂ only by exporting DOM.
- E3 ζ₀ = 1 + linear detachment — no clog, supernatant 0.3 %; established the working set.
- E4 2× field to **60 d** — first steady state (last 10 d: bed +0.14 %, O₂ −0.006 g/m³);
  O₂ 2.96 mg/L consumed, BDOC 89 % (≈27 % total DOC), monotone, surface maximum;
  **H/H₀ plateaus at 1.10 vs Demir's 2.1–4.2 — this set never clogs**; light contrast still
  growing at 60 d (PHO above sand 2.54×) — evidence the response is growth, not trapping.
- E5 µ_PHO ×10/×20 at ζ₀ = 5 — saturates; 40/d ≡ 20/d.
- E6 µ_PHO ×10 on the working set, **run on cosmos** (job 3544074, 4:33 + 1:02) — bed contrast
  1.071, top-2 cm 1.093, supernatant still empty.

**Critique and analysis:** `.claude/CRITIQUE.md` §10–12 (§12 new: buoyant mat export at the top
boundary is a missing mechanism — Nakamoto's mat is lifted by photosynthetic O₂ and leaves via
the overflow, which is why his uncovered filters clog *less* while Campos2002's clog *earlier*;
we are structurally locked into the Campos2002 regime). Two partner-agent picks
(`PICK-A/B-2026-08-26.md`) and the reconciled set in `reports/results-2026-08-26.pdf`.
Literature side-by-side figure `analysis/results/figures/side_by_side_literature_2026-08-26.png`.

**Infrastructure:** cosmos brought current (`src/`, `analysis/` rsynced; corrected preset
verified remotely); `slurm/mu_z01.sbatch` added there.

### Plan for next session

1. **Decide the wet/dry mass convention** (f_dry). It sets whether 1e-3 *is* Campos's field
   influent in our units, moves every µg C/g figure by ~4×, and decides whether the model has a
   single self-consistent operating point. Nothing downstream is settled without it.
2. Submit the **104 d lit/dark pair** at the working set on cosmos (Campos's horizon; contains
   Demir's 55 d) — the only run that can settle the head-loss trajectory and the top-2 cm
   comparison at d97 like-for-like.
3. Decide whether to pursue the clogging behaviour: it needs a higher particulate load (3e-3
   gave H/H₀ 1.29 at 20 d, still short) or the surface-mat mechanisms of CRITIQUE §6/§12.
4. Seasons on the corrected preset (never run) and the OAT at the settled set (still on the
   pre-audit model, at 30 cells, with pathogen removal as its only output).
5. Write the decision file for the linear detachment form (physics argument: shear ∝ v/r
   through a closing pore) — currently in use with no record.

### Open questions / risks

- **Wet vs dry mass** — the state variable is wet (ρ_P 1117 from cell volume × density) while
  RWQM1 stoichiometry is per dry organic matter; yields are ~4× too strong per kg of state
  variable. Fix = ρ_P/f_dry with influents scaled; φ_b, head loss, clogging and detachment are
  invariant, which is also the verification gate.
- **The model cannot clog at realistic loads.** E4 reaches hydraulic steady state at H/H₀ = 1.10
  and stays there; no filter run would ever end.
- **The covered/uncovered contrast is structural, not parametric.** µ_PHO, K_HPO₄, HPO₄_in, δ,
  η_sand and the transfer form have all been swept; the contrast stays at 3–9 % in the bed
  against Campos's 4×. The binding constraint is the ~4 mm lit attaching zone (CRITIQUE §6).
- **Two uncited numbers** in the working set: K_DOM 3e-4 and the ×10 transfer multiplier.
- **`pathogenModel.m:55` MarkerGrowth** still carries the pre-audit half-saturations and drives
  the whole OAT campaign; untouched deliberately, needs its own decision.
- **Operator error:** local `chain_mu_z01_x10_*` files were deleted while assuming they were
  partial when the dark arm had finished; cosmos regenerated them. Verify before deleting.
- **112 untracked files and 13 modified** on `julia-port`, including `simulate.m` and the whole
  preset. A crash or a bad rsync loses two days of work.

## 2026-08-27/28 — autonomous multi-agent session + Manriquez2026 figure recreation

Two Opus 5 agents (model/campaign) with a Fable 5 critic; then a recreation phase. Written
retrospectively by the critic — the per-agent plan files the repo convention asks for were NOT
written before implementation; this entry and the decision file are the record.

**Done** (details: EXPERIMENTS.md E8/E9/P1/P2, `reports/autonomous-2026-08-27.{typ,pdf}`,
`.claude/decisions/2026-08-27-pathogen-model-audit.md`):
- E9: E4 extended to 104 d — steady state terminal, H/H₀ 1.095, Demir unreachable by run length.
  New finding: net N/P mineralisation (effluent NH₄ 8.4×, HPO₄ 6.6× influent), unvalidated.
- E8 load ladder 1×/4× field to 60 d: every axis monotone in load; 4× reaches terminal H/H₀ 1.19
  at O₂ consumption 4.44 — clogging and realistic O₂ not jointly reachable from load.
- pathogenModel audited (MarkerGrowth → audited HET row + HPO₄ term; ζ₀/detachment/flowing-
  inertness now options). P1: 8 pulses on the 60 d filter — the bed passes the pulse
  (L_min 0.07–0.16); PAT export term does not close (5e-5) → log removal defensible only to
  L≈2.9; 100× clogs. OAT campaign staged, not run (`slurm/oat_pulse.sbatch`).
- Figure recreation (`analysis/recreateFigsManriquez2026.m`, figures in
  `analysis/results/figures/recreation/` + README): all results figures in the published design
  on working-set data. Support runs: `fld2x_winter` (90 d, 3 °C, winter light; probeChain gained
  `Temperature`/`LightForm`), `fld2x_cov01` (1 % light, 30 d), 20 d scrapes, day-37 pulses incl.
  the 1e-3× "Pulse 2" (the published BigPulse pairing was 1e-3×, not 100×; 100× clogs at
  day 30.43 and is not reproducible on the working set).
- `analysis/plotBiologicalActivity.m`: reaction-term activity figure, validated against the O₂
  balance (−3.4 %); the filter is overwhelmingly heterotrophic, light adds 2–7 % per term.

**Learned**: the working-set model INVERTS the published seasonal ordering — winter ≥ summer in
the bed at 90 d (rec_biofilm_seasons) — consistent with Bae2023's standing-stock reading and the
probe-chain result the seasonal confound study isolated. O₂ consumption is by design
(decision: the audit made the filter a net O₂ consumer; correct the published caption, not the
model).

**Next**: f_dry decision; diagnose the PAT export non-closure (`simulate.m:921-944`); find a
measured effluent NH₄/HPO₄ to test the mineralisation; submit the OAT once both are settled.
