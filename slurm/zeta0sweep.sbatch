#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J z0sweep
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH --array=0-4
#SBATCH -o slurm/logs/z0_%a.out
#SBATCH -e slurm/logs/z0_%a.err
#
# zeta_0 sweep at N=500, zeta_1 = 0.27, 20 d chained in 10 d legs.
#
# zeta_0 scales the ENTIRE cohesive flux: lambda = zeta_0*phi_b*(1-phi_b), so the
# effective diffusivity D = lambda*Psi'' is linear in it. The three presets
# disagree by four orders of magnitude -- modelPathogen 1e0, pathogenModel 1e2
# (the manuscript value, results.tex:79), modelLund 1e6 -- and results.tex says
# these values were "found through numerical experimentation to obtain sensible
# results", i.e. tuned rather than derived.
#
# The 2026-08-25 deep study concluded that the bed decline is caused by the
# supernatant layer intercepting light and substrate. zeta_0 sets how fast that
# layer spreads, so this sweep tests whether that mechanism is robust to the
# parameter that scales it -- and whether the 1e0 preset behaves differently in
# kind, not just degree.
#
# RISK: the cohesive velocity enters w_v in the CFL bound, so large zeta_0 may
# force very small dt or destabilise. Arms are chained, so a failure still leaves
# leg 1 on disk.
module load matlab/2025b
Z0=(0.1 1.0 10 100 1000)
TG=(0p1 1   10 100 1000)
i=$SLURM_ARRAY_TASK_ID
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
  probeChain(\"z0_${TG[$i]}_n500\", Zeta0=${Z0[$i]}, Zeta1=0.27, NCells=500, Days=20, MaxDt=3e-6, LegDays=10)"
