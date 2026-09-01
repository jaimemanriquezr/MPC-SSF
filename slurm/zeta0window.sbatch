#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J z0win
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH --array=0-2
#SBATCH -o slurm/logs/z0w_%a.out
#SBATCH -e slurm/logs/z0w_%a.err
#
# zeta_0 window search, N=500, zeta_1=0.27, 20 d chained in 10 d legs.
#
# WHY: zeta_0 <= 1 removes the bed peak-and-decline but CLOGS the z=0 cell
# (zeta_0=0.1 at t=3.4 d, zeta_0=1.0 at t=17.0 d, both at CELL=0, the sand
# surface). That is NOT a grid artefact: refining dz 2.5x raises phi_b(0) by only
# 1.10x, so it is converging -- to a value that reaches the 0.99 threshold.
# zeta_0=100 (the manuscript value) does not clog but gives a 1.21x decline, a
# 13 cm supernatant layer, and areal mass above Campos2002's measured band.
#
# zeta_0=10 ran 20 d cleanly at BOTH N=200 and N=500 with only a 1.05x decline,
# so a window plausibly exists between the two failure modes. This brackets it.
#
# MaxDt = 5e-5 (CFL-bound at N=500; identical results to 3e-6 and 26% faster).
module load matlab/2025b
Z0=(3   10   30)
TG=(3   10   30)
i=$SLURM_ARRAY_TASK_ID
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
  probeChain(\"z0w_${TG[$i]}_n500\", Zeta0=${Z0[$i]}, Zeta1=0.27, NCells=500, Days=20, MaxDt=5e-5, LegDays=10)"
