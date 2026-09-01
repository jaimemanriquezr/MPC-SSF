#!/bin/bash -l
# Reproduction of Manriquez2026 fig:pulse-outflow (pathogen.tex:80-98) on the
# WORKING SET: marker concentration at the inflow and the outflow of a 30 d
# mature filter, pulse influx from day 30 to day 32, curves to day 35.
#
#   sbatch slurm/pulse_fig.sbatch          # array 0-3
#
# Submit FROM THE REPO ROOT on cosmos (~/MPC-SSF-manuscript).
#
# FOUR runs, feeding TWO panels (the 1x run appears in both):
#   0  p1x    1x c_ref, nominal rates      -> panel (a) "1x"  AND panel (b) "Experiment 1"
#   1  p100x  100x c_ref, nominal rates    -> panel (a) "100x"
#   2  exp2   1x, d_PAT = 2.0              -> panel (b) "Experiment 2"
#   3  exp3   1x, d_PAT = 2e-6, p_PAT = 2  -> panel (b) "Experiment 3"
#
# The amplitudes are manuscriptExperiments.m:366-367 (Pulse = patNom,
# BigPulse = 100*patNom) and the rate pairs are pathogen.tex's caption verbatim:
# "Experiment 2 increases the inactivation rate to d^PAT_20 = 2 and Experiment 3
# decreases the inactivation and bacterivory rates to d^PAT_20 = 2e-6 and
# p^PAT_20 = 2". No invented rates.
#
# Snapshot: chain_fld2x_lit_leg3.mat = t = 30 d of the E4 chain, which is the
# manuscript's pulse day. Baseline influent PAT is 0 (as the chain was grown);
# the pulse is c_ref = 5.36e-3 kg/m3 (Table B.1, results.tex:21) x the factor.
#
# Wall estimate: 5 simulated days at ~200 s/d on cosmos (measured, job 3549855)
# ~ 17 min per task; the 100x arm ran 1.2x slower than nominal there. 4 h is a
# wide margin.
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 04:00:00
#SBATCH -J pulse_fig
#SBATCH -o slurm/logs/pulse_fig_%A_%a.out
#SBATCH -e slurm/logs/pulse_fig_%A_%a.err
#SBATCH --array=0-3

set -euo pipefail
mkdir -p slurm/logs

i=$SLURM_ARRAY_TASK_ID
case "$i" in
    0) TAG=fig_p1x   ; PF=1   ; DPAT=NaN   ; PPAT=NaN ;;
    1) TAG=fig_p100x ; PF=100 ; DPAT=NaN   ; PPAT=NaN ;;
    2) TAG=fig_exp2  ; PF=1   ; DPAT=2.0   ; PPAT=NaN ;;
    3) TAG=fig_exp3  ; PF=1   ; DPAT=2e-6  ; PPAT=2.0 ;;
    *) echo "unexpected array index $i" >&2; exit 1 ;;
esac

echo "=== pulse_fig tag=$TAG task=$i host=$(hostname) start=$(date -Is) ==="
ls -l analysis/probes/data/chain/chain_fld2x_lit_leg3.mat

module load matlab/2025b

# PulseT0/PulseT1 are offsets from the snapshot time (30 d), so 0..2 is
# "day 30 to 32"; TPost = 5 carries the curves to day 35.
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probePulse(\"$TAG\", Snapshot=\"chain_fld2x_lit_leg3.mat\", \
               PulseFactor=$PF, PulseT0=0.0, PulseT1=2.0, TPost=5.0, NFrames=241, \
               InactivationRate=$DPAT, BacterivoryRate=$PPAT, \
               FlowSurge=1.0, Temperature=19, LightScale=1, \
               NCells=500, Zeta0=1, Zeta1=0.27, DetachForm=\"linear\", \
               KDOM=3e-4, KHPO4=1e-6, TransferScale=10, EtaSand=1500, Delta=5e-3, \
               MaxDt=5e-5, \
               Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])"

echo "=== done $(date -Is) ==="
