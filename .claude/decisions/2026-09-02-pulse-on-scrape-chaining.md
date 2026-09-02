# fig:pat-scraping is re-run by chaining scrape → feed in one process

## The decision

`fig:pat-scraping` is recreated from a new cosmos array (`slurm/patscrape.sbatch`,
job 3563617) rather than kept from April. Each of the four arms (scrape depths
0/4/8/12 cm) runs **one MATLAB process** that:

1. scrapes the end state of leg 3 of `chain_fld2x_lit` (t = 30 d, E4 working set)
   to the arm's depth and runs 20 d of regrowth;
2. saves that end state as a probeChain-format restart snapshot;
3. runs `probePulse` three times on that snapshot — a constant feed at the
   Table B.1 nominal `CRef`, a 10× spike, and a 100× spike, all on `TPost = 3 d`.

To make step 2 possible, `analysis/reproduceManriquez2026.m` gained a
`SnapshotTag` option on `runScrape` (default `""` = the 2026-08-27 behaviour, so
`slurm/repro_scrape.sbatch` is unchanged). It writes `snap` + `rec` via
`stateFromFrame` into `analysis/probes/data/chain/`, with the same `allZero`
assertion `probeChain.m:191` uses.

## Why

**Three feeds, not one.** Jaime asked for two pulse magnitudes (to compare how
challenge magnitude affects removal) plus a constant feed (to show the effect of
feed length). There is also a discrepancy the re-run has to survive either way:
the caption at `pathogen.tex:155` says "a constant feed of marker", while the
revision brief specified a pulse. The 20 d regrowth is the expensive part
(~65 min of the ~2.5 h arm) and is shared across all three feeds, so each extra
feed costs ~25 min. Running all three removes the dependency on settling the
caption before the run.

**One process, not two jobs.** Chaining inside a single task avoids a second
submission that would have to wait on the first array to land. On a deadline
night (co-author draft due 2026-09-02) that latency is the binding constraint,
not CPU.

**Feed order is cheapest-risk-first.** `checkRunFlag` errors on a non-OK flag and
the script runs under `set -euo pipefail`, so an aborting feed kills the task.
The 100× arm is last: if it clogs, the constant-feed and 10× results are already
on disk. This is the failure that lost the ×100 arm of job 3549870 (see the note
at `analysis/probes/probePulse.m`, which now also writes an ABORTED file).

**The pulse runs the same model as the regrowth.** `workingSetSetup` in
`reproduceManriquez2026.m` passes `RespirationForm="reichert"` and `probePulse`
does not. This is a no-op, not a mid-chain model change: every respiration branch
in `src/presets/modelLund.m:134-203` is gated on `PhototrophRespiration > 0`,
which is `0.0` on both sides. The same mismatch already exists in the audited
probeChain → probePulse path, so the new array is consistent with the E-chain.

## Alternatives rejected

- **Keep the April figure with an orange caption note.** It predates the pathogen
  model audit (`.claude/decisions/2026-08-27-pathogen-model-audit.md`), so it
  shows a filter the rest of the revision no longer describes. Retained only as
  the fallback if job 3563617 misses 21:00.
- **Two-stage submission** (scrape array, fetch, then a pulse array). Cleaner
  separation, but adds a queue round-trip and a manual fetch between the stages.
  Rejected on latency, not on correctness.
- **Extend `runScrape` to run the feed itself.** Would duplicate `probePulse`'s
  influent and model construction a second time. The snapshot hand-off keeps
  `probePulse` the single owner of the pulse definition.
- **A `Days = 20` re-run of `repro_scrape.sbatch` alone.** Produces the regrowth
  but no pathogen feed, so it cannot make `fig:pat-scraping`.

## Open

The caption at `pathogen.tex:155` still describes the constant feed only. Which
of the three feeds becomes the figure — and whether the caption gains the two
pulse magnitudes — is Jaime's call once the array lands.
