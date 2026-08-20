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
