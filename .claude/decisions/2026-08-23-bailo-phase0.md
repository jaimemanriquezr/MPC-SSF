# Bailo Phase 0 — is the scheme worth building?

**Nothing of the Bailo scheme is implemented.** Every number below measures the EXISTING
solver. Phase 0 exists to decide whether building Bailo is justified; Bailo results would be
the Phase 1 gates (mass, bounds, monotone energy, second-order accuracy, bounds holding at
100x the explicit dt limit) and none have been run.

Status: **in progress.** 0d complete; 0a–0c (job 3531534) and 0e (job 3531536) running.
Plan: `.claude/plans/2026-08-23-bailo-scheme.md`. Probe: `analysis/probes/probeCflBudget.m`.

## Instrumentation

`simulate.m` gained `RecordCflBudget` (default false), mirroring the existing
`RecordLimitation` pattern. It accumulates exact counters every step — no per-step arrays, so
a million-step run costs nothing — and additionally captures `uCH = xCH(1:n0)`, Solver A's own
φ_b, which the solver normally discards.

**Verified inert when off:** the fixed-step golden reference re-exported **bit-identically**
after each of the two instrumentation changes (`diff -rq` against
`SSF.jl/test/golden/reference`, clean both times).

## 0d — does the negativity guard actually trip? **NO.**

The plan treated bound preservation as justified even at zero speedup, on the grounds that the
guard trips in practice. It does not.

There are two distinct guards, and they are not equivalent for this question:

- `simulate.m:309` **CLOGGED** — `phiBiofilm > CloggingFraction`. This is the upper bound on
  φ_b, i.e. the thing Bailo's bound preservation actually addresses.
- `simulate.m:602,614` **BIOFILM** — an individual *component* concentration `< 0` or `NaN`.
  Bailo bounds the total φ_b and does **not** prevent this.

Evidence:

| source | records | OK | BIOFILM | CLOGGED |
|---|---|---|---|---|
| `analysis/probes/data/*.mat` (66 files) | 64 flags | **64** | 0 | 0 |
| cosmos `slurm/logs` (102 files) | 31 explicit `flag=OK` | 31 | **0** | **0** |

Neither `Unphysical concentration` nor `Clogged!` appears anywhere in 102 cosmos production
logs.

**A stale claim, corrected.** `SSF.jl/test/golden/README.md` says the golden suites avoid the
Lund preset because it "trips a negativity guard immediately". Re-tested directly:
`modelLund()` (κ=1e-6, ζ₀=1e6) at N=100 over 0.5 d returns **`FLAG = OK`**. That claim no
longer reproduces and should not be cited as motivation. It was cited as live evidence in the
`PLANS.md` entry and in conversation on 2026-08-23; both are wrong on this point.

**Consequence.** The "bound preservation is worth it regardless of speed" argument has no
empirical support. Bailo must now be justified either by 0c (a real speed case) or by the
structural component-elimination design, whose own precondition is 0e.

## 0a–0c, 0e — pending

### 0e — the two routes to phi_b differ by O(dt) only

Measured quantity is literally `|u^{n+1} - sum_i phi_i|`: `uCHrec = xCH(1:n0)` is Solver A's
own phi_b (the vector the solver normally discards), against `biofilmFractionLocal` =
`W_enc/rho_L + sum_particles (X_mat + X_enc)/rho_P + sum_liquids S_enc/rho_L` on cells 1..n0.

**Both are sampled at the same time level.** Verified in the loop, not assumed:
`globalBiofilm` is advanced to n+1 at `simulate.m:605`; the frame writes those same
concentrations at `:654`; `uCHrec` is written at `:663` from `:458` in the same iteration.
This mattered — a one-step offset would produce O(dt) drift too, and would have passed the
scaling test below while measuring nothing.

| MaxDt | max abs diff |
|---|---|
| 3e-6 | 2.7281e-05 |
| 7.5e-7 | 6.8254e-06 |

Ratio **3.997** against a dt ratio of 4: first order, so a genuine operator-splitting
inconsistency that converges away, not a structural difference. The two routes are the same
equation, and the component-elimination design is well posed.

Note the probe's original fixed 1e-6 threshold was the wrong test and has been removed: an
O(dt) mismatch is expected at any finite dt, so the scaling is the criterion.

### 0a–0c

Jobs 3531534 (N=100, N=500) and 3531536 (0e, N=100), all κ=1e-6, matured 3 d so the
supernatant carries biofilm. A clean start makes 0c and 0e null tests: at 0.3 d
`phib_sup = 0` and the measured X was 1.000 for reasons that have nothing to do with the
scheme.
