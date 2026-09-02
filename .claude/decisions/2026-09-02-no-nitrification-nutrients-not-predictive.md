# No nitrification: NH4 and HPO4 are nutrient pools, and their effluent values are not predictions

## The decision

Jaime, 2026-09-02. The ecological submodel has **no nitrification and no phosphate sorption**.
NH4 and HPO4 are therefore modelled *purely as nutrients* — pools that reactions draw on and
release — and **their effluent concentrations are not predictions of the model and must be
disregarded**. The manuscript should say so plainly rather than present them as results.

## Why

Measured on the existing chains, influent against effluent (bottom-cell flowing
concentration), summer 19 C at 30 d:

| species | influent | effluent | change |
|---|---|---|---|
| O2 | 9.10e-3 | 6.13e-3 | -33 % |
| IC | 6.23e-3 | 7.25e-3 | +16 % |
| **NH4** | 2.00e-5 | 1.68e-4 | **+739 %** |
| **HPO4** | 5.00e-6 | 3.25e-5 | **+549 %** |
| DOM | 1.00e-3 | 1.14e-4 | -89 % |

The filter **exports** nitrogen and phosphorus. Influent particulate biomass (HET 3.0e-4 plus
PHO 1.0e-3, together larger than the influent DOM at 1.0e-3) attaches, dies and mineralises,
and the released N and P have no sink.

**The literature says the opposite sign**, at comparable influent concentrations.
`documents/Trikkanad2025.pdf` (Trikannad, van der Hoek, Huang, van Halem, *ACS EST Water*
2025, **5**, 6961-6969) reports depth-resolved profiles in full-scale and laboratory SSFs:

- **NH4 is removed completely within the first 45 cm**, and the removal "was accompanied by
  an increase in the NO2- and NO3- concentrations, along with a decrease in pH and DO,
  collectively indicating **nitrification**".
- **PO4 falls from 0.04 mg/L** by 0.02 mg/L in the top 5 cm and a further 0.01 mg/L between
  55 and 90 cm, "reducing PO43- to ultralow levels".
- Their **full-scale** filters had influent NH4 and PO4 **below 0.01 mg/L** — the same order as
  ours (0.02 and 0.005 mg/L) — and still showed "minor but significant decreases". So the
  disagreement is not an artefact of a thin influent.

Two sinks that real filters have and this model does not: **nitrification** (NH4 -> NO2 -> NO3),
and **phosphate sorption** onto sand, which in iron/manganese-coated media removes
ortho-phosphate to below 10 ppb at >99 % efficiency. With neither present, mineralisation is
unopposed and the model is a net mineraliser by construction, not by calibration.

## What this does NOT invalidate

The **internal** nutrient dynamics stand. Phosphate limiting heterotroph growth through the
top 8 cm (E12, E13) is a statement about local Monod terms and the competition between
phototrophs and heterotrophs for a scarce solute; the winter sub-surface bulge was confirmed
by a controlled probe (3565542) that removed the external supply and destroyed the bulge.
That mechanism is unaffected by the absence of a nitrification sink downstream.

What is affected is the **absolute concentration level** those dynamics run at, and every
statement about what leaves the filter.

## Alternatives rejected

- **Add nitrification to the submodel.** The right fix, and a genuine model extension: two
  more reactions, two more components (NO2, NO3), new stoichiometry and new half-saturations.
  Out of scope for this revision, which is already carrying an audited kinetics change.
- **Tune the death rates down until the export disappears.** Fitting a symptom. The export is
  structural — a source with no sink — and suppressing it would corrupt the biomass dynamics
  to fix the solute chemistry.
- **Say nothing and keep the outflow figures as they are.** Rejected because
  `fig:1d-outflow-liquids` and `fig:2d-plots` both include NH4 and HPO4 panels, so silence
  would leave the paper asserting effluent values the depth-resolved literature contradicts.

## Consequences for the manuscript

- `fig:1d-outflow-liquids` (Figure 8) and `fig:2d-plots` (Figure 9) each show five species
  including NH4 and HPO4. Those panels need an explicit caveat, or the two species dropped.
- Any discussion sentence that reads the effluent nutrient concentrations as a result must go.
- The model's standing as a *biofilm* model is unaffected; this is a limitation of the
  ecological submodel's solute chemistry, and stating it is stronger than leaving a reviewer
  to find it.
