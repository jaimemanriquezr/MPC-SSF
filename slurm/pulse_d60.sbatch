#!/bin/bash -l
# Figure 11 (fig:pulse-outflow) extended from t = 37 d to t = 60 d.
#
#   sbatch slurm/pulse_d60.sbatch          # array 0-3, one per panel arm
#
# The existing panels stop at t = 37 d (TPost = 7) because of an earlier wall-time
# limit, not because the physics stops there. This runs the SAME four arms for
# TPost = 30, i.e. t = 30 -> 60 d, so the figure can be drawn at 40, 50 or 60 d
# from ONE run per arm rather than three.
#
# The four arms are the published panel set, parameters read off the rec of the
# existing pulse_fig_*_d37 files rather than assumed:
#   p1x    PulseFactor 1      preset rates (inact 0.02, bacterivory 8.0)
#   p1em3  PulseFactor 1e-3   preset rates
#   exp2   PulseFactor 1      InactivationRate 2      (high inactivation)
#   exp3   PulseFactor 1      InactivationRate 2e-6, BacterivoryRate 2 (low both)
# All four feed on [0, 2) d after the snapshot, i.e. days 30-32, at CRef.
#
# RESOLUTION. NFrames = 1441 keeps the existing 30-minute sampling across 30 d
# (the d37 files used 337 frames over 7 d, the same interval). ~314 MB per arm,
# 1.2 GB for the set; cosmos has 135 TB free.
#
# Wall: the d37 runs took ~22 min for 7 d = 3.2 min per simulated day, so ~96 min
# for 30 d. 4 h covers an arm that clogs and drives the adaptive step down --
# note pulse_fig_p100x ABORTED at t = 30.42 d, so a long marker run CAN fail;
# probePulse writes an ABORTED file in that case rather than losing the trajectory.
#
# Output: analysis/probes/data/pulse/pulse_fig_{p1x,p1em3,exp2,exp3}_d60.mat
# recreateFigsManriquez2026's pulseFile prefers _d37 over the plain stem; it will
# need a _d60 entry added to pick these up.
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH -J pulsed60
#SBATCH -o slurm/logs/pulse_d60_%A_%a.out
#SBATCH -e slurm/logs/pulse_d60_%A_%a.err
#SBATCH --array=0-3

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) STEM=p1x;   EXTRA="PulseFactor=1"                                             ;;
    1) STEM=p1em3; EXTRA="PulseFactor=1e-3"                                          ;;
    2) STEM=exp2;  EXTRA="PulseFactor=1, InactivationRate=2"                         ;;
    3) STEM=exp3;  EXTRA="PulseFactor=1, InactivationRate=2e-6, BacterivoryRate=2"   ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== pulse_d60 stem=$STEM extra=[$EXTRA] host=$(hostname) start=$(date -Is) ==="
ls -l analysis/probes/data/chain/chain_fld2x_lit_leg3.mat

module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probePulse(\"fig_${STEM}_d60\", Snapshot=\"chain_fld2x_lit_leg3.mat\", NCells=500, \
        PulseT0=0.0, PulseT1=2.0, TPost=30.0, NFrames=1441, $EXTRA)"

echo "=== done $(date -Is) ==="
