#!/bin/bash -l
# Submit the 20 d / kappa=1e-6 / N=500 pair. Job B reruns the same script after
# A finishes (afterany, not afterok: if A dies on wall time having completed
# summer, B still picks up winter via the E3_longterm/done.txt sentinel).
set -euo pipefail
cd "$HOME/MPC-SSF-manuscript"
mkdir -p slurm/logs
A=$(sbatch --parsable slurm/manuscript_k1e6_20d.sbatch)
B=$(sbatch --parsable --dependency=afterany:"$A" slurm/manuscript_k1e6_20d.sbatch)
echo "submitted A=$A (summer+winter)  B=$B (resumes winter if A ran out of time)"
squeue -u "$USER" -o '%.10i %.9P %.8j %.2t %.10M %.10l %R'
