# Biomass sand-attachment factor: can surface concentration reproduce Schijven's scrape?

## The hypothesis

E13 found the model cannot reproduce Schijven2013's 0.6-1.6 log scraping loss at ANY marker
sand factor: dL is 0.02 log at best, because removal is spread over 2-25 cm and the top 2 cm
carries only 18 % of it. The model does not concentrate removal at the surface; real filters
evidently do.

The marker's own attachment cannot fix this -- E13 swept it and the depth distribution did
not move, because the marker does not build the filter. What builds the filter is HET and
PHO attachment, and both currently attach to bare sand at full strength
(SandAttachmentFactor = 1).

**Hypothesis.** Lowering the BIOMASS sand factor makes heterotrophs and phototrophs attach
preferentially where biofilm already exists rather than to clean grains. That is a positive
feedback: biofilm accretes on biofilm. It should thicken the Schmutzdecke, thin the deep bed,
and so make a 2 cm scrape cost more -- moving dL towards Schijven's 0.6-1.6 log.

Done means: 20 d chains grown at several biomass sand factors, the depth distribution of
biofilm measured, and dL re-measured intact vs scraped for each.

## Sweep

`SandBiomass` in {0.1, 0.3, 1.0}. 1.0 is the current model. Marker stays at s = 0 (the
published choice) so the two effects are not confounded; E13 already characterised the
marker axis independently.

## Steps

1. `src/presets/modelLund.m` — add `options.SandBiomass (1,1) {mustBeNumeric} = 1.0`, applied
   as the `SandAttachmentFactor` of HET and PHO only. POM is unaffected in practice (its
   AttachmentSand is 0, so its factor is inert). Default 1.0 leaves every existing caller
   byte-identical. Verify: `modelLund()` gives HET/PHO factor 1; `SandBiomass=0.3` gives 0.3.
2. `src/presets/pathogenModel.m` — forward `SandBiomass` (NaN = preset), as for Attenuation.
3. `analysis/probes/probeChain.m` and `probePulse.m` — same NaN-forwarded option. probePulse
   MUST carry it too: it rebuilds the model, and a challenge run on a chain grown at a
   different factor would be a mongrel.
4. `slurm/sandbio_chain.sbatch` — 3 arms, 20 d chains, tags `sb01`, `sb03`, `sb10`.
5. `slurm/sandbio_scrape.sbatch` — 6 arms, {intact, scraped} x 3, submitted with
   `--dependency=afterok` on stage 1. Each scraped arm builds its OWN uniquely named
   snapshot, so the arms cannot race.

## Files touched

`src/presets/modelLund.m`, `src/presets/pathogenModel.m`, `analysis/probes/probeChain.m`,
`analysis/probes/probePulse.m`, two sbatch scripts.

## What would falsify the hypothesis

If the biofilm depth profile is unchanged by SandBiomass, the feedback does not operate and
the surface concentration must come from somewhere else -- straining, which this model does
not represent at all, being the obvious candidate. That is a publishable negative too:
it would say a continuum attachment model cannot produce a Schmutzdecke effect of the
measured size without an explicit straining term.
