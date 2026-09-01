#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J deepstd
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH --array=0-9
#SBATCH -o slurm/logs/deep_%a.out
#SBATCH -e slurm/logs/deep_%a.err
#
# Deep study, 2026-08-25. Balanced across three questions:
#
#   NUMERICAL   tasks 0,1  -- timestep cap at 20 d (1e-6 reference, 5e-5 loose)
#               task  9    -- N=1000 refinement confirmation
#               task  8    -- Dirichlet vs Neumann at z = 0
#   MECHANISM   tasks 2-5  -- detachment sweep 0 / 0.5 / 1 / 2. THE test for the
#               bed peak-and-decline: loss ~ phi_b*sqrt(q/(1-phi_b)) is
#               super-linear in phi_b and strongest at the bed top (Leg A). If the
#               peak survives DetachScale = 0, the decline is instead the
#               supernatant layer starving the bed of light and substrate (Leg B).
#   VALIDATION  every task saves effluent time series for all nine species, so
#               removal can be compared with Bellamy1985 / Bae2023 / Elemo2024.
#
# All runs N=500 (or 1000), chained in 10 d legs: each leg writes its own
# chain_<tag>_leg<k>.mat, so a killed or over-long job still yields partial
# results, and resubmitting resumes from the last completed leg.
#
# Chaining verified by probeRestart: no snapshot field is zeroed, and a chained
# run matches a continuous one to ~1e-4 relative (residue = adaptive stepper
# re-ramping at the leg boundary, suppressed by seeding AdaptiveInitialDt).
module load matlab/2025b
TAG=(  dt1e6   dt5e5   det0    det0p5  det1    det2    z001    z015    bcdir      n1000  )
Z1=(   0.27    0.27    0.27    0.27    0.27    0.27    0.01    0.15    0.27       0.27   )
DET=(  1.0     1.0     0.0     0.5     1.0     2.0     1.0     1.0     1.0        1.0    )
NC=(   500     500     500     500     500     500     500     500     500        1000   )
DAY=(  20      20      20      20      20      20      20      20      20         10     )
DT=(   1e-6    5e-5    3e-6    3e-6    3e-6    3e-6    3e-6    3e-6    3e-6       3e-6   )
BC=(   neumann neumann neumann neumann neumann neumann neumann neumann dirichlet  neumann)
i=$SLURM_ARRAY_TASK_ID
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
  probeChain(\"${TAG[$i]}\", Zeta1=${Z1[$i]}, DetachScale=${DET[$i]}, NCells=${NC[$i]}, \
             Days=${DAY[$i]}, MaxDt=${DT[$i]}, CohesionBC=\"${BC[$i]}\", LegDays=10)"
