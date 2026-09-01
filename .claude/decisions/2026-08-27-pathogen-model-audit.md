# `pathogenModel` audited against the manuscript and the corrected Lund set (2026-08-27)

Requested: "Audit `pathogenModel.m` the way `modelLund` was audited." Companion to
`.claude/decisions/2026-08-25-half-saturations-corrected.md`, which corrected the Lund
half-saturations and explicitly deferred the marker reactions:

> **Correcting `pathogenModel.m:55` MarkerGrowth (O₂ 3e-3, NH₄ 4e-3, DOM 2e-4 — the old
> HET set) at the same time.** NOT changed: it controls PAT growth and therefore every
> number in the OAT/Sobol campaign behind `sensitivity.tex`. Needs its own decision.

This is that decision. Sources read: `manuscripts/AWR-SSF/results.tex:128-165`
(`tab:eco-parameters`, the table the docstring calls "Table B.4"),
`manuscripts/AWR-SSF/pathogen.tex:20-60` (the marker reaction vector),
`src/@State/simulate.m:1115-1137` (how Monod terms are combined),
`.claude/decisions/2026-08-25-half-saturations-corrected.md`,
`.claude/decisions/2026-08-25-kinetics-wolf-reichert.md`.

## Finding 0 — the Lund base is inherited, not duplicated

`src/presets/pathogenModel.m:36` calls `modelLund(...)` and mutates the returned arrays.
No Lund constant is written twice, so the 2026-08-25 audit already reached the pathogen
preset. **Verified at runtime, not by reading**: `analysis/testPathogen.m` part (a)
compares every component field (`Name`, `Density`, `Dispersivity`, `TransportRate`,
`Attenuation`, `AttachmentSand`, `AttachmentMatrix`, `SandAttachmentFactor`) and every
Lund reaction field (`NominalRate`, `TemperatureCorrectionFactor`, `IsLightDependent`,
`MinimumLightFactor`, `OptimalLightFactor`, `LightInhibition`, `EfficiencyBiofilm`, and
every half-saturation and stoichiometric key) one by one, prints the comparison, and
asserts that exactly **three** differences exist, all deliberate:

| difference | modelLund | pathogenModel | why |
|---|---|---|---|
| PAT `TransportRate` | 5.47 /d | 218.8 /d (×40) | `srun_pathogen` convention, carried from the Julia campaign |
| PAT `SandAttachmentFactor` | 1 | 0 (`SandPathogen`) | `pathogen.tex:17` sets `b^{att,PAT}_f = 0` |
| every Lund `EfficiencyFlowing` | 1 | 0 | Lund reactions inert in the flowing suspension (see Finding 4) |

Printed output (`matlab -batch "addpath('analysis'); testPathogen"`):

```
Heterotroph growth  K_O2    0.0002    0.0002      Phototroph growth  K_IC    0.0012  0.0012
Heterotroph growth  K_DOM   0.004     0.004       Phototroph growth  K_NH4   2e-05   2e-05
Heterotroph growth  K_NH4   1e-06     1e-06       Phototroph growth  K_HPO4  2e-05   2e-05
Heterotroph growth  K_HPO4  2e-05     2e-05       Phototroph growth  rate    2       2
Heterotroph growth  rate    2         2           Phototroph death   rate    0.276   0.276
Hydrolysis          K_POM/HET 0.1     0.1         Hydrolysis         rate    3       3
```

## Finding 1 — `tab:eco-parameters` gives the marker NO half-saturations at all

The table lists exactly three marker constants — µ_PAT 2.00e-1, d_PAT 2.00e-2,
p_PAT 8.0 /d (`results.tex:140-142`), all with reference "---" (own calibration). There
is **no** K^*_PAT row, **no** K_pred row, and **no** θ^PAT row, although `pathogen.tex:55-59`
declares θ^PAT_{20,µ}, θ^PAT_{20,d} and θ_{20,pred} as model constants. The implementation
therefore had to borrow, and what it borrowed was the **pre-audit heterotroph-growth set**.

Since `pathogen.tex:26` and the stoichiometry `[PAT +1, O2 −1.2317, DOM −1.5873]` make r7
literally the heterotroph growth row with PAT in place of HET, the borrowing is right —
it just borrowed the superseded values. **Decision: r7's Monod set is the audited
heterotroph-growth set of 2026-08-25, and stays tied to it** (the test asserts equality
with `modelLund`'s heterotroph row key by key, so a future Lund correction cannot leave
the marker behind again — which is exactly how this drift happened).

## Finding 2 — HPO₄ was missing from the marker Monod set

`manuscripts/AWR-SSF/pathogen.tex:33`:

> `M^PAT_* := min_{* in {O2, NH4, HPO4, DOM}} monod(*)`

Four terms, combined with **min**. The code had three (`O2, NH4, DOM`). The solver does
use min — `src/@State/simulate.m:1129`, `[monod(:, i), j] = min(mon_term, [], 2)` — so
this was a straight omission of a term the manuscript's own equation requires, not a
modelling choice. **HPO₄ added at the audited K^HPO4_HET = 2.0e-5.**

## Constants before / after

| constant | before | after | source |
|---|---|---|---|
| MarkerGrowth K_O₂ | 3.0e-3 | **2.0e-4** | Reichert2001 RWQM1 K_O2,H = 0.2 g/m³ (= Wolf = ASM); via the audited HET row |
| MarkerGrowth K_NH₄ | 4.0e-3 | **1.0e-6** | Wolf2007 K_S,H,NH₃ ≈ 0 ("N never limits"); 1e-6 is the depletion-protective stand-in agreed on 2026-08-25 |
| MarkerGrowth K_HPO₄ | *absent* | **2.0e-5** | Reichert2001 K_HPO4,H = 0.02 g P/m³; required by `pathogen.tex:33` |
| MarkerGrowth K_DOM | 2.0e-4 | **4.0e-3** | Wolf2007 K_S,H,SS |
| MarkerGrowth µ_PAT | 0.2 /d | **0.2 /d** (unchanged) | `results.tex:140`, own calibration |
| MarkerGrowth θ | 1.047 | **1.047** (unchanged) | not in the table; equals the published θ^HET_{20,µ} = θ^PHO_{20,µ} = 1.047 (`results.tex:143-144`) |
| Inactivation d_PAT | 0.02 /d | **0.02 /d** (unchanged) | `results.tex:141` |
| Inactivation θ | 1.08 | **1.08** (unchanged) | not in the table; equals the published θ^PHO_{20,d} = 1.080 |
| Bacterivory p_PAT | 8.0 /d | **8.0 /d** (unchanged) | `results.tex:142` |
| Bacterivory θ | 1.08 | **1.08** (unchanged) | not in the table |
| Bacterivory K_pred | 2.0e-3 | **2.0e-3** (unchanged) | **the table gives no value**; the `thesis_model` number is kept, as the old docstring already said |
| Bacterivory `EfficiencyFlowing` | 1.0e-3 | **1.0e-3** (unchanged) | see Finding 3 — the manuscript says 1.0e-2 |
| PAT `TransportRate` | ×40 | ×40 (unchanged) | `results.tex:118` gives b^{tr,PAT}_P = 2.18e2 = 40 × 5.47 ✓ **confirmed by the table** |
| PAT `SandAttachmentFactor` | 0 | 0 (unchanged) | `results.tex:112` b^{att,PAT}_sand = 0 ✓ **confirmed** |
| Cohesion ζ₀ | 1e2, hardcoded | **1e2, now the `Zeta0` option default** | `results.tex:79` is the manuscript value; the E1–E7 working set is 1 |
| Detachment | `1.4e-5*sqrt(|v|/18)`, hardcoded | **same, now `DetachForm="campaign"` default**, plus `"sqrt"` (0.14·√) and `"linear"` (0.14·|v|/18) | the two named forms are `probeChain.m`'s `detachmentFunction` verbatim |

Every default is the previous value, so all 40 existing call sites of `pathogenModel`
are unchanged apart from the four MarkerGrowth constants.

## Finding 3 — flowing-phase bacterivory efficiency disagrees with the manuscript

`pathogen.tex:60` states `E^PAT_{f,p} = 0.01`. The code's `WaterFactor` default is
`1e-3`. **NOT changed**: it is the value the completed Julia OAT/Sobol campaign behind
`sensitivity.tex` ran at, and changing the default would silently move every number in
it. Recorded in the `pathogenModel` docstring as a caveat with the fix (`WaterFactor=1e-2`).
Someone has to decide which of the two the manuscript will claim.

## Finding 4 — the Lund reactions' flowing-phase inertness (raised by the reviewer)

Objection received: running a pulse on an E4 snapshot with `EfficiencyFlowing = 0` on the
Lund reactions means "the mature filter is not the same model that grew it".

**The premise is false, and it is checkable in one line.** `analysis/probes/probeChain.m:77`
is `mp = pathogenModel(...)`. probeChain builds *its* model from this preset too, so the
E1–E7 chain snapshots — `chain_fld2x_{lit,dark}_leg{1..6}` included — were themselves grown
with every Lund reaction inert in the flowing suspension. Keeping the inertness during the
pulse is what preserves the host model; switching it off would be the mid-chain change.

**Decision: keep it on by default, and expose it.** New option `LundFlowingInert`
(default `true`). `false` restores `modelLund`'s own flowing efficiencies, so the contrast
is a one-argument run for anyone who wants to measure it. `testPathogen` asserts both
directions. The quantified 3 d lit comparison the reviewer asked for was **not run** — the
comparison the reviewer proposed measures the difference between the snapshot's own physics
and a different model, which is a question about probeChain's design, not about whether the
pulse matches its host. It is a legitimate separate experiment and is not claimed here.

## Finding 5 — new open finding: the PAT mass budget does not close to solver tolerance

Found while writing the mass-balance test, not part of the audit brief.

With **every reaction rate set to zero** and a constant influent, the ε-weighted budget
`dM = supply − export` (the convention of `analysis/probes/probeMassClosure.m`) closes to

| species | residual / supply | exported? |
|---|---|---|
| HET | 2.94e-11 | no (export 3.3e-10) |
| PHO | 2.94e-11 | no (export 1.1e-9) |
| **PAT** | **4.96e-5** | **yes (export 1.17e-2, 61 % of supply)** |

Frame refinement 101 → 501 → 2001 gives −6.27e-5 → 4.44e-5 → 4.86e-5, converging to a
nonzero value, so it is not trapezoid error. The distinguishing feature is that PAT is the
only particle that actually leaves the column. `probeMassClosure.m` ran with influent
PAT = 0 and therefore **never exercised its own outflow term** — its 1.1e-9 result is a
statement about retained species only, and `.claude/CRITIQUE.md` §10's "§2.1 is CLOSED"
should be read with that restriction.

Likely cause: the budget's export term is `q·c_out`, while `src/@State/simulate.m:936,944`
exports `porosityBoundaries(end)·volumeAvgVelocity(end)·globalFlowing(end)`. The two agree
only if `ε·v_avg = q` exactly at the outflow face, which biofilm growth and osmosis need
not preserve. **Not diagnosed further here.** `analysis/testPathogen.m` asserts absolute
closure at 1e-4 (not 1e-6) and states why in its header; the conservative-tracer property
is instead tested as *invariance* — the PAT residual is −5.5659e-6 with the Lund biology
running and −5.5654e-6 with every reaction off (1.04e-4 relative), which is the real
statement that no Lund reaction touches PAT.

## Alternatives rejected

- **Leaving MarkerGrowth on the pre-audit set for campaign continuity.** That is what the
  2026-08-25 deferral did, and the cost is that the marker's O₂ half-saturation (3 mg/L)
  sat an order of magnitude above the heterotrophs' (0.2 mg/L) for two days in a model
  where r7 is *defined* as the heterotroph growth row. Continuity is preserved the honest
  way instead: the OAT campaign must be re-run (it must be re-run anyway — its nominals
  are the pre-audit Lund rates; see `analysis/PARAMETERS.md`).
- **Adding HPO₄ as an option defaulting to off.** Rejected: `pathogen.tex:33` is not
  optional, and a switch would preserve the divergence it is meant to remove.
- **Changing `WaterFactor` to the manuscript's 1e-2.** Deferred, see Finding 3.
- **Changing the marker θ values to something cited.** Nothing to cite — the table has no
  PAT θ rows. Left alone and recorded as such.

## Consequences to carry

- **Every PAT result produced before today used K_O₂ = 3e-3, K_NH₄ = 4e-3, K_DOM = 2e-4
  and no HPO₄ term.** That is the whole Julia OAT/Sobol campaign behind `sensitivity.tex`.
- With the min-Monod and the added HPO₄ term, **r7 is identically zero whenever influent
  HPO₄ is zero** — which is the case for the manuscript influent (`results.tex` Table B.1,
  HPO₄ = 0) and for `logOatCampaign.m:35`. The `marker_growth` OAT parameter therefore has
  exactly zero sensitivity at that influent. This is not new behaviour introduced by
  choice: heterotroph growth already had a K_HPO₄ term and was already dead there. It is
  an argument for running the campaign at the field influent (HPO₄ = 5e-6), which is what
  `slurm/oat_pulse.sbatch` does.
- `tab:eco-parameters` needs new rows for K^*_PAT and K_pred, or a sentence saying r7 uses
  the heterotroph constants. At present the paper's Monod set for PAT is unreproducible
  from the paper.
- SSF.jl mirrors these presets; its `pathogen_model` needs the same four constants and the
  goldens re-anchoring. Not done here (separate repo).
