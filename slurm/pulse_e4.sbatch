#!/bin/bash -l
# PAT challenge pulses on the E4 60 d mature filter (chain_fld2x_{lit,dark}_leg6),
# audited pathogenModel, working set (zeta0 = 1, zeta1 = 0.27, LINEAR detachment
# x1, transfer x10, K_DOM 3e-4, K_HPO4 1e-6, eta_sand 1500, delta 5 mm, N = 500).
#
#   sbatch slurm/pulse_e4.sbatch          # array 0-7
#
# Submit FROM THE REPO ROOT on cosmos (~/MPC-SSF-manuscript); every addpath is
# relative to it. `#!/bin/bash -l` is required or `module` is not found.
#
# Four operating conditions x {lit, dark}:
#   nom    nominal 10x pulse   (peak 5.36e-2 kg/m3 on [0.1, 0.3) d)
#   big    BigPulse 100x       (peak 5.36e-1)
#   surge  10x pulse + hydraulic surge x1.4 on the inflow velocity
#   cold   10x pulse at 3 C -- THETA-CORRECTED KINETICS ONLY. The biofilm in the
#          snapshot was grown at 19 C, so this is a cold shock of a summer
#          filter, not a winter filter. Do not report it as a winter run.
#
# PREREQUISITE: chain_fld2x_{lit,dark}_leg6.mat must be in
# analysis/probes/data/chain/ on cosmos (they are -- job 3549014 resumes from
# them). probePulse errors out with probePulse:noSnapshot if not.
#
# Wall estimate: 3 simulated days at ~92 s/d locally (measured, 0.35 d smoke run
# = 32 s at N = 500), cosmos ~1.8x slower (E6) => ~9 min per task. The BigPulse
# and 3 C arms are stiffer; 4 h is a wide margin.
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 04:00:00
#SBATCH -J pulse_e4
#SBATCH -o slurm/logs/pulse_e4_%A_%a.out
#SBATCH -e slurm/logs/pulse_e4_%A_%a.err
#SBATCH --array=0-7

set -euo pipefail
mkdir -p slurm/logs

# index:      0        1         2        3         4          5           6         7
ARM=(lit      dark     lit       dark     lit       dark       lit         dark)
COND=(nom     nom      big       big      surge     surge      cold        cold)

i=$SLURM_ARRAY_TASK_ID
A=${ARM[$i]}
C=${COND[$i]}
case "$A" in
    lit)  LIGHT=1 ;;
    dark) LIGHT=0 ;;
esac
case "$C" in
    nom)   PF=10  ; SURGE=1.0 ; TEMP=19 ;;
    big)   PF=100 ; SURGE=1.0 ; TEMP=19 ;;
    surge) PF=10  ; SURGE=1.4 ; TEMP=19 ;;
    cold)  PF=10  ; SURGE=1.0 ; TEMP=3  ;;
esac
TAG="e4_${A}_${C}"

echo "=== pulse_e4 tag=$TAG task=$i host=$(hostname) start=$(date -Is) ==="
ls -l analysis/probes/data/chain/chain_fld2x_${A}_leg6.mat

module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probePulse(\"$TAG\", Snapshot=\"chain_fld2x_${A}_leg6.mat\", \
               PulseFactor=$PF, FlowSurge=$SURGE, Temperature=$TEMP, \
               LightScale=$LIGHT, TPost=3.0, NFrames=145, \
               NCells=500, Zeta0=1, Zeta1=0.27, DetachForm=\"linear\", \
               KDOM=3e-4, KHPO4=1e-6, TransferScale=10, EtaSand=1500, Delta=5e-3, \
               MaxDt=5e-5, \
               Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])"

echo "=== done $(date -Is) ==="
