# Julia consolidated into SSF.jl; MPC-SSF is MATLAB-only (2026-08-22)

## The decision

The Julia implementation now lives solely in `code/1d/SSF.jl`. `MPC-SSF/julia`
is retired, the package is renamed `MPCSSF` -> `SSF`, and the Julia batch scripts
move with the code. MPC-SSF is MATLAB-only from here.

Approved by Jaime 2026-08-22 ("I want the duplication gone").

## Why

The Julia sources were duplicated: `MPC-SSF/julia` and the split-out package,
kept in sync by hand through the `julia-split` subtree branch. The cost was
concrete and had already been paid twice:

- The two trees had **diverged**. `MPC-SSF/julia` carried the cardinal (CTMI)
  temperature response (`src/ecological/Reaction.jl`, commit 89558ec) that the
  standalone package did not: 283 tests vs 276. A rename applied to one tree
  would have conflicted across ~29 files on the next merge.
- Nine of thirteen `slurm/*.sbatch` scripts were Julia jobs sitting in what is
  now a MATLAB repo, `cd`-ing into a directory that no longer exists.

Keeping one Julia home removes the sync obligation entirely.

## Evidence (measured, not inferred)

- Trees compared before deletion: 363 vs 362 files, the difference being a
  gitignored `Manifest.toml`. Six files differed; only `Reaction.jl` and
  `runtests.jl` were substantive (the CTMI work).
- After merging `julia-split` into SSF.jl: **283/283 tests pass**, matching what
  `MPC-SSF/julia` ran, with all four golden reference suites bit-identical
  (`RESULT: MATCH - Julia reproduces the MATLAB pathogen reference`).
- `julia/Manifest.toml` (a stray nested lockfile) was the only tracked file
  unique to `MPC-SSF/julia`; SSF.jl had deliberately removed it in `fce9968`.
- Slurm account defect: all nine Julia scripts carried `-A lu2025-7-124`. Both
  that and `lu2026-2-100` are valid for `jrman` on `lu48` (checked with
  `sacctmgr show assoc user=jrman`), so the wrong account billed silently.

## Alternatives rejected

- **Keep both trees, apply changes twice.** This is the status quo that produced
  the 283/276 divergence. It fails whenever one side is edited alone.
- **Delete `MPC-SSF/julia` without merging first.** Would have discarded the
  CTMI implementation from the working tree (recoverable from `julia-port`
  history, but absent from the surviving package) and four files of Jaime's
  uncommitted path fixes (not recoverable at all).
- **Merge `matlab-claude` and expect `julia/` to disappear.** It does not:
  `julia/test/golden/compare_pathogen.m` exists only on that branch, so a merge
  treats it as an addition and keeps it. It was relocated to SSF.jl instead.
- **Regenerate the package UUID on rename.** Julia resolves by UUID, not name,
  so keeping `72a1329e-...` leaves existing environments valid. A fresh UUID is
  the right call only if the package is ever registered, where `SSF` will also
  face a name-collision review.
- **Point the moved slurm scripts at a new absolute cosmos path.** Absolute
  paths are exactly what broke twice this month. They now rely on Slurm starting
  in the submit directory, matching the MATLAB scripts.

## Consequences

- CLAUDE.md's cross-implementation source of truth moves to
  `code/1d/SSF.jl/test/golden/`.
- The cosmos trees (`~/MPC-SSF*`) are rsync copies, not git repos, and all four
  still contain `julia/`. They are stale but harmless; nothing there depends on
  the rename. Redeploy when a Julia run is next needed, letting cosmos resolve
  its own Manifest (it has Julia <= 1.10.4; the local one is 1.12.6).
- SSF.jl's GitHub remote was renamed to `SSF.jl`; the local URL was updated.
