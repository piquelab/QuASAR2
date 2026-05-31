#!/bin/bash
#SBATCH -J quasar2_sim
#SBATCH -q primary
#SBATCH --mem=160G
#SBATCH --time=4-00:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -o logs/quasar2_sim_%j.out
#SBATCH -e logs/quasar2_sim_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=go7535@wayne.edu

set -euo pipefail

mkdir -p logs

module unload gnu7
module load gnu9 R

export LD_LIBRARY_PATH=/wsu/home/groups/piquelab/apps/el7/misc/lib:/wsu/el7/groups/piquelab/R/4.5.1/lib64/R/lib:${LD_LIBRARY_PATH:-}

cd /rs/rs_grp_scaipgenetic/QuASAR2/QuASAR2ex

echo "Started at: $(date)"
echo "Running on: $(hostname)"

Rscript quasar2_simulations_power_analyses.R

echo "Finished at: $(date)"