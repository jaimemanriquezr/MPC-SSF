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
