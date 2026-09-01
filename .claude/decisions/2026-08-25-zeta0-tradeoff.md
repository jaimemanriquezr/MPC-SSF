# ζ₀ is a tradeoff dial, not a fix; ζ₀ = 5 adopted as the working value (2026-08-25)

## The decision

Working parameter set for all further manuscript-scenario runs: **ζ₀ = 5, ζ₁ = 0.27,
κ = 1e-6, δ = 5 mm, N = 500, ImplicitDispersion = true, AdaptiveMaxDt = 5e-5**, shin
scheme, Neumann cohesion BC. Chosen as the best point on a monotone tradeoff curve,
NOT as a value that removes the bed peak-and-decline.

## Why (measured)

Runs: `analysis/probes/data/chain/chain_z0w_{3,10,30}_n500_leg{1,2}.mat` (cosmos,
20 d), `chain_z0_100_n500_leg{1,2}.mat` (cosmos, 20 d), `chain_z0_1_n500_leg1.mat`
(cosmos, 10 d; leg 2 aborted CLOGGED at t ≈ 17.0 d, no file), `chain_z0w_5_n500_local_leg{1..6}.mat`
(local, 40 d in legs 3/2/5/10/10/10 d), `chain_z0w_5_n1000_local_leg1.mat` (local, 10 d).
All ζ₁ = 0.27, N = 500, MaxDt = 5e-5, detachment 0.14·√(|v|/18). Head loss from
`analysis/headlossKozenyCarman.py` → `analysis/results/headloss_zeta0_2026-08-25.txt`.

| ζ₀ | t_end | bed peak | at | bed(20 d) | decline(20 d) | sup % (20 d) | φ_b(0) (20 d) | H/H₀ 10 d → end | t(H/H₀=3) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 10 (clogged 17) | 0.06084 | 10.0 | — | — | 5.1 (10 d) | 0.7997 (10 d) | 7.07 | 4.5 d |
| 3 | 20 | 0.06651 | 19.8 | 0.06650 | 1.00× | 22.7 | 0.8128 | 6.08 → 5.41 | 4.8 d |
| **5** | **40** | 0.06213 | 16.8 | 0.06137 | 1.01× (1.12× at 40 d) | 30.0 (41.5 at 40 d) | 0.7453 (0.845 at 40 d) | 5.68 → 3.57 | 4.7 d |
| 10 | 20 | 0.05815 | 13.8 | 0.05485 | 1.06× | 38.9 | 0.6582 | 5.24 → 3.41 | 4.8 d |
| 30 | 20 | 0.05443 | 12.8 | 0.04700 | 1.16× | 49.1 | 0.5441 | 4.75 → 2.57 | 5.0 d |
| 100 (paper) | 20 | 0.05243 | 11.8 | 0.04198 | 1.25× | 55.5 | 0.4550 | 4.47 → 2.18 | 5.0 d |

- Every column is monotone in ζ₀: lower ζ₀ delays the turnover and reduces the decline
  but loads the surface cell. No interior optimum exists.
- Every arm declines; "no decline" at ζ₀ = 3/5 in the 20 d table is a window artefact
  (ζ₀ = 5 loses 1.12× by 40 d).
- ζ₀ = 5 does **not** clog: φ_b(0) saturates at 0.845 (rate +0.0008/d at 37–40 d).
  Total biomass reaches steady state by 40 d (+0.25 % over the last 5 d); bed still
  declining slowly (−0.0002/d), supernatant 41.5 %.
- Convergence at ζ₀ = 5: N = 1000 vs 500 at 10 d, bed +0.37 %, profile L2 0.25 %,
  φ_b(0) 0.5853 vs 0.5877. N = 500 is converged in this regime too.
- Head loss (Kozeny–Carman, uniform bed): H/H₀ crosses 3× at 4.5–5.0 d for **every**
  ζ₀ and then falls. ζ₀ does not control operational clogging; the influent load does.
- Wall (N = 500, per 10 d): 631 s (10–20 d), 1175 s (20–30 d), 1550 s (30–40 d) — cost
  per simulated day rises with φ_b. N = 1000: 2349 s for 0–10 d.

## Alternatives rejected

- ζ₀ ≤ 1 (monotone bed growth, Campos-like): clogs at N = 500 (0.1 at 3.4 d, 1 at
  17.0 d) while surviving at N = 200 — under-resolution, not physics. Unusable.
- ζ₀ = 3: same decline behaviour as 5 within the window, 7 points more surface load
  (0.813 vs 0.745 at 20 d) and no 40 d evidence that it saturates.
- ζ₀ = 100 (manuscript): earliest turnover, 1.25× decline, 55 % of biomass above sand.

## Caveats

- The ζ₀ = 5 arm is six chained legs; before any number goes in the paper re-run it as
  one continuous job (`probeChain` chaining is gated to 4.36e-5 by `probeRestart`).
- All ζ₁ evidence is at ζ₀ = 100; the ζ₁ sweep must be repeated at ζ₀ = 5.
- See `.claude/CRITIQUE.md` §3, §10 for what this parameter does NOT fix.
