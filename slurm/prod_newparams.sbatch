#!/bin/bash -l
# PRODUCTION chains under the new parameters, 90 d, both influent phosphate values.
#
#   sbatch slurm/prod_newparams.sbatch          # array 0-3
#
# New presets in force (commit b1fcd03 + 9c637bf): mu_HET = 4.8 (Tenore2021
# Table 1 mu_max,2), nu_P = 52 (Tenore k_tot = 210 on a wet-mass basis). Every
# existing chain was grown at 2.0 and 0.094, so nothing in the manuscript has yet
# been produced with these.
#
# BOTH influent phosphate values are run because the choice is open:
#   p0  HPO4_in = 0        the PUBLISHED tab:incoming-concentrations value
#   p5  HPO4_in = 5.0e-6   the 2026-08-25 field set, what the working set uses
# E12 showed HPO4 = 0 costs ~5 % of areal biomass at 60 d and destroys the winter
# sub-surface bulge, so the two are not interchangeable. And the filter is a net
# phosphate SOURCE -- effluent runs 6.5x influent at 19 C by mineralisation of
# influent biomass -- so HPO4_in = 0 does not starve it.
#
# Summer and winter only; the dark/covered arm can follow once these are read.
#
# Wall: summer ~3.2 min per simulated day at N = 500 (measured, d37 pulses), so
# 90 d ~ 4.8 h. Winter is FASTER, ~1.3 min/d measured on 3565542 (slow kinetics
# permit larger adaptive steps), so ~2 h. 12 h covers either with margin. Legs are
# written individually, so a wall-clock kill loses at most one leg and the job can
# be resubmitted to continue.
#
# Output: analysis/probes/data/chain/chain_prod_{lit,win}_{p0,p5}_leg{1..9}.mat
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -J prodnew
#SBATCH -o slurm/logs/prod_newparams_%A_%a.out
#SBATCH -e slurm/logs/prod_newparams_%A_%a.err
#SBATCH --array=0-3

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) TAG=prod_lit_p0; HPO4=0.0;    TEMP=19; LFORM=chain  ;;
    1) TAG=prod_lit_p5; HPO4=5.0e-6; TEMP=19; LFORM=chain  ;;
    2) TAG=prod_win_p0; HPO4=0.0;    TEMP=3;  LFORM=winter ;;
    3) TAG=prod_win_p5; HPO4=5.0e-6; TEMP=3;  LFORM=winter ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== prod_newparams tag=$TAG HPO4_in=$HPO4 T=$TEMP host=$(hostname) start=$(date -Is) ==="
echo "--- presets in force ---"
grep -n "attenuationParticle = \|NominalRate=4.80" src/presets/modelLund.m || true

module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probeChain(\"$TAG\", Zeta0=1, Zeta1=0.27, NCells=500, Days=90, MaxDt=5e-5, LegDays=10, \
               FramesPerLeg=241, LightScale=1, LightForm=\"$LFORM\", Temperature=$TEMP, \
               KDOM=3e-4, TransferScale=10, KHPO4=1e-6, DetachForm=\"linear\", \
               Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, $HPO4, 1.0e-3])"

echo "=== done $(date -Is) ==="
