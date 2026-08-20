# Winter/summer biomass inversion: accepted, resolved by seasonal influent (2026-08-20)

## The decision

The winter > summer standing-biomass inversion under season-constant influent is
ACCEPTED as a valid model output; the E1/E2 protocol adopts a SEASONAL influent
(winter PHO_in = 1e-4 kg/m3 = 20% of the 5e-4 summer field value, both within
Campos2006b Fig. 1(a) seasonality). Approved by Jaime 2026-08-20.

## Why (measured, not inferred) — probeBudget.m, data/budget_*.mat

Closure budget dM/dt = q(c_in - c_out) + growth - death - endog over 10 d,
field influent, full corrected structure (rates, half-sats, endogenous
respiration, source-faithful death):

- Growth is NOT minor: ~1.7x supply in both seasons (the earlier trap-and-decay
  reading was wrong on this point). Export ~0: all influent algae are captured.
- Per-capita, summer grows 1.8x faster (0.31 vs 0.17 /d) — but per-capita LOSSES
  are 3.5x seasonal (0.34 vs 0.10 /d) because theta_loss = 1.08 > theta_growth
  = 1.047 and self-shading/nutrients buffer growth's temperature response.
  Net per-capita: summer -0.03 /d, winter +0.075 /d. Cold suppresses decay
  harder than growth: that is the whole mechanism.
- Winter stock is linear in supply (1.69 kg per kg/d*10d, three-point sweep
  5e-4/1e-4/2e-5). Flip threshold: winter influent < 49% of summer. Real
  seasonality (5-25%) is far below it — no tuning involved.
- Field echo: winter bed-deposit accumulation is observed (Nakamoto via
  Campos2006b p.888 — winter clogging despite lower influent SS).

## Alternatives rejected

- Fixing via rates/structure: probed exhaustively (PG-excess, endogenous, death
  0.4 -> 0.09) — the ratio never moved; it lives in the theta asymmetry, which is
  source-faithful.
- Filamentous surface-growth mode (Campos A_in): not needed for the flip; stays
  future work for mechanistic summer-mat formation.

## New audit entry

Our theta_growth,PHO = 1.047 is Campos's BACTERIA value; their algae theta_kga =
1.066 (same a/b swap as their Table 3 symbol listing). Softens but cannot reorder
theta_growth vs theta_loss; recorded for the Table 2 corrections.

## Addendum: winter light-curve lever REFUTED (2026-08-20, late; probeLight.m)

Sweep of the winter light amplitude (peaks 0.6/0.3/0.15, 10 d, decided structure,
theta_growth,PHO = 1.066): a 4x light cut changes winter growth by 6% and stock by
5% (0.0593 -> 0.0566 kg/m2; W/S 1.81 -> 1.73). Two causes: (1) Steele concavity —
l(0.15) = 0.35 vs l(0.6) = 0.90, so light cuts are sub-linear in the factor;
(2) growth is RECYCLED-PHOSPHORUS-limited, not light-limited: influent HPO4 = 0,
so growth tracks the P released by death (+0.0037) and endogenous respiration
(+0.01) — biomass-proportional, light-independent fluxes. The loop
biomass -> losses -> P recycling -> growth makes the seasonal ratio pure theta
arithmetic; under fixed influent, no growth-side parameter (mu, theta, light) can
flip it. Only identified fixed-influent lever: route decay-released P into POM
(hydrolysis-gated) instead of dissolved HPO4 — a structural change. This
strengthens the adopted seasonal-influent resolution.

## Addendum: the inversion is mostly the respiration term (probeNoResp.m)

Clean pair with endogenous respiration removed (death 0.09 O2-neutral the only
biological loss): W/S collapses 1.81 -> 1.10 (summer PHO 0.076, winter 0.083).
Respiration drives the flip via (1) ~75% of the theta-sensitive loss pool and
(2) its P recycling (+0.01/unit vs death's +0.0037) fueling winter's
light-independent growth. Residual 1.10 = trap-and-decay on death alone; the
PG-excess-era inversion (3.6x) had a different cause (dark biomass creation).
Explicit tension: the term that delivers the O2 deficit (without it effluent
returns to 9.2-9.4, supersaturated) amplifies the winter inversion under fixed
influent — deficit and summer-dominance cannot coexist with season-constant
water. The adopted seasonal influent resolves both simultaneously.
