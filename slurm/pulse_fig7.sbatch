#!/bin/bash -l
# fig:pulse-outflow RECREATION runs (2026-08-28) — supersedes pulse_fig.sbatch's
# data for the recreated panels. Two corrections against the 08-27 set:
#   * the published fig:pat-low-inact pairs the 1x pulse with a 1e-3x amplitude
#     (Pulse 2 inflow plateau ~5e-6 g/L in PATPulse_HighLow.pdf), NOT the 100x
#     BigPulse of manuscriptExperiments E5 — so task 1 is PF=1e-3 (no clogging);
#   * the published time axis runs to day 37, so TPost=7 (was 5).
#
#   0  p1x_d37    1x   c_ref, nominal rates   -> Pulse 1 in both panels
#   1  p1em3_d37  1e-3 c_ref, nominal rates   -> Pulse 2, panel (a)
#   2  exp2_d37   1x, d_PAT = 2.0             -> Experiment 2, panel (b)
#   3  exp3_d37   1x, d_PAT = 2e-6, p_PAT = 2 -> Experiment 3, panel (b)
#
# From chain_fld2x_lit_leg3.mat (t = 30 d), pulse day 30->32, curves to day 37.
# ~20 min per task at cosmos speed. Submit FROM THE REPO ROOT.
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 3:00:00
#SBATCH -J pulseFig7
#SBATCH -o slurm/logs/pulse_fig7_%A_%a.out
#SBATCH -e slurm/logs/pulse_fig7_%A_%a.err
#SBATCH --array=0-3

set -euo pipefail
mkdir -p slurm/logs

i=$SLURM_ARRAY_TASK_ID
case "$i" in
    0) TAG=fig_p1x_d37   ; PF=1     ; DPAT=NaN  ; PPAT=NaN ;;
    1) TAG=fig_p1em3_d37 ; PF=1e-3  ; DPAT=NaN  ; PPAT=NaN ;;
    2) TAG=fig_exp2_d37  ; PF=1     ; DPAT=2.0  ; PPAT=NaN ;;
    3) TAG=fig_exp3_d37  ; PF=1     ; DPAT=2e-6 ; PPAT=2.0 ;;
    *) echo "unexpected array index $i" >&2; exit 1 ;;
esac

echo "=== pulse_fig7 tag=$TAG task=$i host=$(hostname) start=$(date -Is) ==="
module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probePulse(\"$TAG\", Snapshot=\"chain_fld2x_lit_leg3.mat\", \
               PulseFactor=$PF, PulseT0=0.0, PulseT1=2.0, TPost=7.0, NFrames=337, \
               InactivationRate=$DPAT, BacterivoryRate=$PPAT, \
               FlowSurge=1.0, Temperature=19, LightScale=1, \
               NCells=500, Zeta0=1, Zeta1=0.27, DetachForm=\"linear\", \
               KDOM=3e-4, KHPO4=1e-6, TransferScale=10, EtaSand=1500, Delta=5e-3, \
               MaxDt=5e-5, \
               Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])"

echo "=== done $(date -Is) ==="
