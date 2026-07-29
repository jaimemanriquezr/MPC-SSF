# o2pen is a censored metric; report NaN and add o2dep alongside it

## The decision

`qois()` in `morris_screening.jl` returns `o2pen = NaN` whenever the oxygen
threshold is not crossed strictly inside the domain, instead of the depth of the
last cell. A new QoI `o2dep` — the fraction of oxygen consumed across the domain
— is added to `QOI_NAMES` to carry the signal `o2pen` was meant to carry.

## Why

The Morris screen reported $mu^* = 0$ *exactly* for every parameter except
`velocity` on `o2pen`. That is not physical insensitivity, it is saturation.

`o2pen` was

```julia
o2pen_i = findlast(>(0.01 * maximum(o2) + 1e-30), o2)
o2pen   = o2pen_i === nothing ? 0.0 : z[o2pen_i]
```

`findlast` returns the last index at which oxygen exceeds 1% of its in-domain
maximum. If the filter never approaches anoxia, that index **is the last cell**,
for every parameter set, and the metric returns the domain bottom as a constant.

Measured at `ncells = 30`, mature biofilm, `pulse`, 1 day:

| quantity | value |
|---|---|
| max O2 | 1.093e-2 |
| min O2 | 6.601e-3 |
| min / max | **0.604** |
| cells at or below the 1% threshold | **0 of 62** |
| `findlast` index | 62 of 62 → z = 1.0164 m = domain bottom |

Oxygen falls only to 60% of its maximum. The 1% threshold is never approached, so
the metric cannot resolve anything. For contrast the HET penetration metric `pen`,
which uses the same 1%-of-surface construction, crosses at cell 44 of 62 and
behaves normally — the flaw is in this threshold against this profile, not in the
pattern.

Returning the domain bottom is the dangerous failure: it is a plausible depth, it
is dimensionally correct, and it reads as "oxygen reaches the bottom of the
filter" when the truthful statement is "this metric never resolved". `NaN`
propagates through the elementary-effect means and makes the column visibly
unusable instead of quietly constant.

## Why add o2dep rather than only fix o2pen

`NaN` alone would leave the screen with no oxygen QoI at all, and oxygen
consumption is a real quantity worth ranking. `o2dep = 1 - min(O2)/max(O2)` is
continuous, bounded, and defined whether or not the filter goes anoxic, so it
cannot saturate the way a threshold-crossing depth does.

It varies usefully. Across the `velocity` range at otherwise nominal settings:

| velocity (m/d) | `o2dep` | `o2pen` |
|---|---|---|
| 2.16 | 0.4249 | censored |
| 7.2 | 0.3962 | censored |
| 21.6 | 0.5000 | censored |

Note the direction: consumption *rises* with velocity despite the shorter
residence time, because faster flow drives biomass deeper — `pen` moves from
0.393 m to 0.885 m over the same range. That is a real result the old metric
could not have produced.

## Alternatives rejected

**Lower the threshold from 1% to something reachable, e.g. 70% of maximum.** Then
the depth resolves, but the number means "where oxygen has dropped 30%", which is
not a penetration depth and would be reported under a name implying it is. The
censoring is a property of the profile, not of the constant.

**Leave `o2pen` returning the bottom and document the caveat.** Rejected: a caveat
in a plan file does not travel with a CSV column, and the column looks like data.

**Drop the oxygen QoI entirely.** Simplest, but discards a measurable quantity for
which a well-behaved metric exists.

## Consequence

`QOI_NAMES` grows from 9 entries to 10. It is used only inside
`morris_screening.jl`, so nothing downstream needed changing. The existing
`results/morris/mu_sigma.csv` is pre-fix and stale for the separate cohesion
reason, and is being regenerated anyway.
