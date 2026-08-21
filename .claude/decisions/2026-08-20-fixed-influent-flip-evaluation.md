# Fixed-influent flip evaluation — G-probe and closure findings (2026-08-20/21)

Plan: ~/.claude/plans/unified-juggling-widget.md (adversarial, Opus critic, 2 rounds).
Probes: analysis/probes/{probePInventory,probePClosure}.m, data/{pinv,pclo}_*.mat.

## G-probe (P inventory/flux decomposition, critic-designed)

- tau_P = 0.05-0.07 d BOTH seasons -> approach A (P-routing through POM) is DEAD:
  its precondition (tau_winter >> 10 d >> tau_summer) fails by 3 orders of magnitude.
- P CLOSURE VIOLATED (~45% of P input created): HET death releases 0.0209 P per
  unit vs 0.0141 consumed by growth (N: 0.0653 vs 0.0248 — 2.6x). HET death is
  90% of all P release; the growth loop runs substantially on phantom P. PHO
  death under-releases (0.0037 vs 0.01), destroying P. NEW audit-family finding:
  the Lund death rows do not conserve nutrients against their own growth rows.
- Production/retention separated: summer produces MORE gross PHO (0.33 vs 0.23
  kg/m2 per 30 d) but winter RETAINS (stock W/S 2.03 at theta_growth 1.066).

## H (nutrient-closed death rows): bug fix confirmed, flip worsened

Dissolved release at death set equal to biomass content (HET N 0.0248 P 0.0141;
PHO N 0.06 P 0.01; POM carries none). Closure defect drops 10x (0.0040->0.0004
summer). BUT W/S rises 2.03 -> 2.51: the phantom P was a theta-1.08-scaled
(summer-biased) subsidy — removing it costs summer more (gross PHO -26% summer,
-13% winter). ADOPT as a correctness fix regardless of the flip; it joins the
Table 2 correction list. It does not rescue the fixed-influent route.

## Cardinal (CTMI) machinery: implemented, tested, both ports

Reaction TemperatureResponse="cardinal" + CardinalTemperatures=[Tmin Topt Tmax],
mu(20 C) = NominalRate anchor. MATLAB testCardinal.m passes; Julia 283/283 incl.
new testset; goldens untouched. E ensemble BLOCKED on the pre-registered
Bernard & Remond 2012 triplet set (paper paywalled/bot-blocked; requested from
Jaime into documents/).

## Standing

Fixed-influent candidates: A dead, C dead (round 1), F likely-null (no run),
B conditional-on-E, H adopted-but-worsens, E pending triplets. If E fails on
honest triplets, the fixed-influent route closes with pre-registered evidence
and the seasonal influent stands alone.

## Addendum 2026-08-21: theta_growth sweep pre-empts E — evaluation closing

Jaime asked why the TCF itself was never swept; diagnostic sweep run (cosmos
array 3528283, probeTheta.m, data/theta_*.mat): theta_growth,PHO 1.066/1.09/1.12
gives W/S = 1.81/1.76/1.72 — a 2.3x suppression of winter mu moves realized
winter growth by 8%. Same self-compensation as the light sweep: growth is
P-SUPPLY-limited, so HPO4 accumulates and the Monod climbs when mu drops.

Consequence: in the mu_winter -> 0 limit, winter stock -> supply/losses ~ 0.036
vs summer ~ 0.033: W/S ~ 1.1 (the trap-and-decay floor), still >= 1. The earlier
prediction that a cardinal triplet flips to W/S ~ 0.28 was WRONG (mishandled the
supply feedback); no growth-side intervention of ANY strength can flip under
fixed influent. E is pre-empted regardless of triplet — the Bernard & Remond
blocker is moot. Confirmation run (winter, theta=2.0, mu-factor 8e-6) submitted
as cosmos job 3528328; on COMPLETED+consistent, the fixed-influent route is
CLOSED and the seasonal influent stands as the unique resolution (B's only
escape — a nutrient-independent mat — was already shown to require an external
nutrient boundary condition, i.e., seasonal boundary again).

## Addendum 2026-08-21b: asymptote run exposes the budget-residual method

Winter theta=2.0 (mu-factor 8e-6, growth literally off; cosmos 3528329):
PHO stock 0.0539, but the budget still attributes 0.0489 "growth" — a
mu-INDEPENDENT residual. The residual method (probeBudget onward) is therefore
measuring a systematic bookkeeping imbalance (~0.05 winter / presumably ~0.065
summer), not photosynthetic growth. Supply convention checked in simulate.m
(flux = porosity * (q/porosity) * c_in = q*c_in — correct); suspicion moves to
phase-vs-bulk concentration conventions in the mass integrals. Decisive check
submitted (cosmos 3528356): zero-biology run where dM must equal supply-export
exactly, plus summer theta=2.0 for the true growth attribution.

WHAT STANDS regardless (direct stock measurements, no residuals involved):
- W/S = 1.81 / 1.76 / 1.72 / 1.65 at theta_growth 1.066 / 1.09 / 1.12 / 2.0.
- The mu -> 0 winter asymptote is MEASURED: W/S floor = 0.0539/0.0326 = 1.65,
  far above 1 — the fixed-influent closure conclusion (no growth-side
  intervention can flip) holds and is STRONGER than the arithmetic predicted.
WHAT IS QUARANTINED pending the mass check: all "growth residual" magnitudes
and the per-capita growth/loss numbers derived from them (incl. the
"P self-compensation" interpretation and the earlier "winter net-growing
+0.075/d" mechanism statement — winter stock is now shown to be essentially
growth-independent, i.e., pure trapping/decay, which REinstates the original
H1 reading).

## Addendum 2026-08-21c: mass check — TRANSPORT CREATES PARTICLE MASS

Zero-biology run (cosmos 3528356_0, probeMassCheck, data/mck_winter_masscheck.mat):
all reaction rates 0, initial state all-zero (State.m verified), export ~0,
q = 7.2 (SandFilter default verified). Hard bound: stock <= cumulative supply.
Measured at 10 d: PHO 0.0536 vs supply 0.0360 (+49%); HET 0.4507 vs 0.1930
(+133%). Particle mass is created in transport. Non-proportionality (stock
ratio 8.41 vs influent ratio 5.36 under identical transport parameters and a
shared velocity field) implies a concentration-dependent inconsistency —
suspect phase-vs-bulk concentration conventions in the flowing/attachment
coupling (simulate.m fluxFlowing/attachment terms), not a constant factor.

Consequences:
- The mu-independent "growth residual" (~0.05) is explained: it was this
  transport excess, not growth. probeBudget-family residuals are dead.
- W/S seasonal ratios compare like-for-like (identical influent both seasons)
  and likely remain qualitatively valid; ABSOLUTE stocks are suspect.
- Scope unknown: the refactor may share this with the legacy solver (published
  results) and/or the Julia mirror (anchored to MATLAB, so it would reproduce
  the same behavior). NEEDS a dedicated investigation: influent-off decay
  test, single-cell analytic comparison, grid-refinement, legacy comparison.
- Summer theta=2.0 asymptote (3528356_1): stock 0.0286 — with growth off and
  the transport excess present, summer/winter absolute numbers are both
  inflated; W/S floor at mu~0 recomputes to 0.0539/0.0286 = 1.88.
