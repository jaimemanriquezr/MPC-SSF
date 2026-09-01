#!/bin/bash -l
# Working-set OAT, STARTUP scenario (clean IC; ripening is the disturbance -- no snapshot).
# on the E1-E7 working set, anchored on the E4 60 d mature filter
# (chain_fld2x_lit_leg6), pulse scenario.
#
#   sbatch slurm/oat_startup.sbatch          # array 0-28, one PARAMETER per task
#
# Submit FROM THE REPO ROOT on cosmos (~/MPC-SSF-manuscript).
#
# NOT YET RUN as of 2026-08-27. Only the ParamFilter smoke tests were executed
# (locally): N=30/TMature 0.5/TPost 0.6 on `marker_growth`, and N=500 from the
# snapshot with TPost 0.35 on `bacterivory`. Both completed end to end.
#
# ONE TASK PER PARAMETER, not per (parameter, sign). logOatCampaign runs both
# the x2 and the x1/2 arm inside one call and needs the unperturbed baseline in
# the same process to form D_min, I_min and the asymmetry, so a task is
# 3 simulations: baseline, plus, minus. Splitting by sign would require a
# baseline in each half anyway and would double the baseline count from 29 to 58.
#
# WALL ESTIMATE. 3 simulated days at N=500 costs ~92 s/d locally (measured:
# probePulse 0.35 d smoke run = 32 s) and cosmos is ~1.8x slower (E6), so
# ~8.5 min per simulation, ~26 min per task, and 29 tasks in parallel finish in
# ~30 min wall / ~12.5 h aggregate CPU. The 4 h limit covers a task that clogs
# and drives the adaptive step down. Perturbed arms that abort are recorded as
# clog_driver = 1 in measures.csv, not as failures.
#
# Results: analysis/results/oat/log_oat_startup_<param>/ (L0.csv, curves_*.csv,
# measures.csv). Concatenate the 29 measures.csv into one table afterwards.
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 04:00:00
#SBATCH -J oat_startup
#SBATCH -o slurm/logs/oat_startup_%A_%a.out
#SBATCH -e slurm/logs/oat_startup_%A_%a.err
#SBATCH --array=0-28

set -euo pipefail
mkdir -p slurm/logs

# The 29-parameter design of logOatCampaign/campaignParams (velocity excluded,
# a fixed operating condition -- decision 2026-08-18). Order must match that
# function; a name not in the design makes the task fail loudly
# (logOatCampaign:unknownParam), it no longer produces an empty measures.csv.
PARAMS=(temperature influent_PAT dispersivity transport_P attach_sand sand_pathogen \
        mu_HET mu_PHO d_HET d_PHO hydrolysis theta_growth theta_death \
        K_O2_HET K_DOM_HET K_HPO4_HET K_O2_PAT K_pred \
        marker_growth inactivation bacterivory \
        zeta_0 kappa zeta_1 detach_scale \
        light_att_water light_att_sand attenuation_P beta_porosity)

P=${PARAMS[$SLURM_ARRAY_TASK_ID]}

echo "=== oat_startup param=$P task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is) ==="
ls -l analysis/probes/data/chain/chain_fld2x_lit_leg6.mat

module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    logOatCampaign(\"startup\", Preset=\"workingset\", \
        NCells=500, \
        TPost=3.0, MaxDt=5e-5, NFrames=145, \
        ParamFilter=\"$P\", OutTag=\"_$P\")"

echo "=== done $(date -Is) ==="
