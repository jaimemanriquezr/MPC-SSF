#!/bin/bash -l
# Can we reproduce Schijven2013's scraping loss? Sand-attachment sweep, intact vs scraped.
#
#   sbatch slurm/sandfactor_scrape.sbatch        # array 0-7
#
# Schijven2013 (D1 -> D2) scraped the Schmutzdecke off a ripened SSF and measured
# 0.6 log10 lower removal of MS2 and 1.6 log10 lower removal of E. coli WR1.
#
# Our marker has SandAttachmentFactor s = 0: it cannot attach to bare sand, only to
# biofilm. Schijven's own model has a sand-attachment BASELINE plus a Schmutzdecke
# term, so s = 0 deletes the baseline he fits. This sweeps s and measures the
# scraping loss Delta L = L(intact) - L(scraped) against his 0.6-1.6 log.
#
# PREDICTION. At s = 0 all removal is biofilm-mediated, so scraping should cost the
# MOST. As s rises, bare sand carries more of the removal and scraping should cost
# LESS. Whichever s brackets 0.6-1.6 log is the value consistent with Schijven.
#
#   s = 0     current model, marker cannot stick to clean sand
#   s = 0.06  Schijven2013 T4 sticking efficiency (geo-mean); the OAT nominal
#   s = 0.3   intermediate
#   s = 1.0   sand as sticky as biofilm -- the Diehl2025 formulation, no asymmetry
#
# Base filter: chain_fld2x_lit_leg2, the 20 d working-set filter. Scraped arm uses
# scrape2cm_from20d.mat, the SAME filter with the top 2 cm removed and NO regrowth
# (reproduceManriquez2026 Days=0), so the pair differs only by the scrape --
# Schijven's D1/D2 design. 2 cm is the depth fig:pat-scraping uses and carries 18 %
# of removal activity (measured 2026-09-02).
#
# CHALLENGE: constant feed at the Table B.1 nominal CRef over the whole 3 d window,
# the closest analogue to Schijven's continuous seeding. Read rec.Lmin per arm.
#
# CAVEAT for any conclusion: Schijven measured MS2 and E. coli; our marker is
# neither, and its inactivation/bacterivory rates are Manriquez Table B.4, not
# fitted to his organisms. A match in Delta L is evidence about the ATTACHMENT
# STRUCTURE, not a validation of the pathogen submodel.
#
# Wall: ~30 min per arm (3 d at N = 500 on cosmos); 8 arms in parallel. 2 h covers
# an arm that clogs.
#
# Output: analysis/probes/data/pulse/pulse_sf_{intact,scraped}_s{0,006,03,10}.mat
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH -J sfscrape
#SBATCH -o slurm/logs/sandfactor_scrape_%A_%a.out
#SBATCH -e slurm/logs/sandfactor_scrape_%A_%a.err
#SBATCH --array=0-7

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) SNAP=chain_fld2x_lit_leg2.mat;  S=0;    TAG=intact_s0    ;;
    1) SNAP=chain_fld2x_lit_leg2.mat;  S=0.06; TAG=intact_s006  ;;
    2) SNAP=chain_fld2x_lit_leg2.mat;  S=0.3;  TAG=intact_s03   ;;
    3) SNAP=chain_fld2x_lit_leg2.mat;  S=1.0;  TAG=intact_s10   ;;
    4) SNAP=scrape2cm_from20d.mat;     S=0;    TAG=scraped_s0   ;;
    5) SNAP=scrape2cm_from20d.mat;     S=0.06; TAG=scraped_s006 ;;
    6) SNAP=scrape2cm_from20d.mat;     S=0.3;  TAG=scraped_s03  ;;
    7) SNAP=scrape2cm_from20d.mat;     S=1.0;  TAG=scraped_s10  ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== sandfactor_scrape tag=$TAG snap=$SNAP s=$S host=$(hostname) start=$(date -Is) ==="
ls -l analysis/probes/data/chain/$SNAP

module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probePulse(\"sf_$TAG\", Snapshot=\"$SNAP\", NCells=500, SandPathogen=$S, \
        PulseFactor=1, PulseT0=0.0, PulseT1=3.0, TPost=3.0)"

echo "=== done $(date -Is) ==="
