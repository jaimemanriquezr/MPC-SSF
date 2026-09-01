# Film-form liquid transfer implemented, then rejected as default: it puts the biofilm maximum inside the bed (2026-08-26)

## The decision

`SolverOptions.TransferForm = "film"` (b = D/L_f², L_f = εφ_b/a_s on grains, mat slab
above the sand; exact relaxation per step) stays in the code but is **not** the default and
is not used for any reported run. Runs use the constant form (`Liquid.TransportRate`,
600/d) with at most a ×10 multiplier. Jaime, 2026-08-26: "No reasonable simulation should
have φ_b peak deep into the sand bed."

## Why

With the phases equilibrated (film form, or constant ×100 and above) the φ_b and
heterotroph maxima move to 2–5 cm below the sand surface
(`chain_fa_film_lit`, `chain_xf_100_lit`, `chain_xf_1000_lit`); at ×1 and ×10 the profile
is monotone from z = 0 (`chain_kd_3e4_lit`, `chain_xf_10_lit`). Campos2002 and Demir2017
put the biomass maximum in the uppermost 1.5–2 cm, decreasing with depth. Mechanism: once
the biofilm's pore water equilibrates with the flow, DOM produced inside the biofilm by
hydrolysis at the surface leaks out, is advected down and consumed over the DOM-uptake
length; slow exchange keeps the hydrolysate where it was made.

The derivation (transfer = a_s·D/L_f, i.e. a diffusion time L_f²/D; 600/d implies a
540 µm film) is not wrong; what is wrong is treating internally produced DOM like O₂, a
small molecule that exchanges freely. Real biofilms retain hydrolysate in the EPS matrix.

## Alternatives

- Constant ×10 (adopted for the operating point): monotone profile, O₂ in the Elemo window,
  BDOC removal 74–86 % (Campos-consistent as a fraction of total DOC), 2.5× cheaper than ×100.
  Uncited multiplier; the D/L_f² argument (L_f ≈ 170 µm) is its provenance.
- Film form + retention of hydrolysis products in the enclosed phase (a retained fraction
  on the hydrolysis source, or a species-specific exchange with DOM slow and O₂/IC fast):
  the physically complete fix; future work, needs both ports.

## Consequences

- `fa_film_*` and the first launch of `delta2cm_*`/`mu_x*` (film form) are diagnostics
  only; the δ/μ arms were killed and relaunched on constant ×10.
- Diagnostics still running to pin the mechanism: `fa_film_lit_n500` (grid dependence of
  the peak position), `fa_film_dom0_lit` (hydrolysate-only feeding).

## Addendum (same day): mechanism attribution revised

`fa_film_lit_n500`: peak at 2.2 cm (2.0 cm at N=200) — grid-independent, a physical length.
`fa_film_dom0_lit` (DOM_in = 0): peak returns to z = 0. So the in-bed maximum is fed by
INFLUENT DOM under free exchange, not by leaked hydrolysate as first stated; the ~2 cm is
the bed's DOM-uptake length. Why the surface cell itself loses the maximum (highest
detachment velocity? trapped-phototroph load?) is not established. The rejection stands
on the profile alone (Campos2002/Demir2017: maximum in the top 1.5–2 cm).
