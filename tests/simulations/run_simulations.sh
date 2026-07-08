#!/bin/bash
# run_simulations.sh
# Master orchestration script for QuASAR2 benchmarking simulations.
#
#   Step 1 — Generate simulated datasets (array: one task per seed)
#   Step 2 — Fit each method in parallel (four arrays, each depends on Step 1)
#   Step 3 — Aggregate results and produce plots (depends on all Step 2 arrays)
#
# USAGE:
#   sbatch run_simulations.sh [--seeds 10] [--N_lo 60] [--N_hi 300] [--M 100] [--delta 0.10]
#
# OUTPUT LAYOUT:
#   results/
#     data/   raw simulated data RDS files
#     CR/     QuASAR2CR result RDS files
#     GLM/    QuASAR_GLM result RDS files
#     Q2/     QuASAR2 result RDS files
#     LM/     LM result RDS files
#   logs/     SLURM stdout/stderr logs
#
# SBATCH directives for the orchestration job itself 
#SBATCH -q primary
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o logs/quasar2_sim_%j.out
#SBATCH -e logs/quasar2_sim_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=go7535@wayne.edu

set -euo pipefail

mkdir -p logs

echo "Started at: $(date)"
echo "Running on: $(hostname)"

# ============================================================
# ENVIRONMENT
# ============================================================
module unload gnu7
module load gnu9 R

# ============================================================
# DEFAULT PARAMETERS  (override via command-line flags)
# ============================================================
N_SEEDS=10
N_LO=60
N_HI=300
M=100
DELTA=0.10

# Paths
SCRIPT_DIR="/rs/rs_grp_scaipgenetic/QuASAR2/tests/simulations"
RESULTS_DIR="/rs/rs_grp_scaipgenetic/QuASAR2/tests/simulations/results"
LOG_DIR="/rs/rs_grp_scaipgenetic/QuASAR2/tests/simulations/logs"
RSCRIPT="Rscript"

# SLURM resource settings for child jobs
QUEUE="primary"
MEM_DATA="8G"
MEM_METHOD="16G"
MEM_AGG="8G"
TIME_DATA="01:00:00"
TIME_METHOD="04:00:00"
TIME_AGG="01:00:00"
CPUS=1
MAIL_USER="go7535@wayne.edu"

# ============================================================
# PARSE COMMAND-LINE ARGS
# ============================================================
while [[ $# -gt 0 ]]; do
  case "$1" in
    --seeds)  N_SEEDS="$2";  shift 2 ;;
    --N_lo)   N_LO="$2";     shift 2 ;;
    --N_hi)   N_HI="$2";     shift 2 ;;
    --M)      M="$2";        shift 2 ;;
    --delta)  DELTA="$2";    shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

DELTA=$(printf "%g" "${DELTA}")
M=$(printf "%g" "${M}")
N_LO=$(printf "%g" "${N_LO}")
N_HI=$(printf "%g" "${N_HI}")

COND_TAG="N${N_LO}-${N_HI}_M${M}_d${DELTA}"

echo "============================================================"
echo "QuASAR2 simulation run"
echo "  Seeds:   1 – ${N_SEEDS}"
echo "  N range: [${N_LO}, ${N_HI}]"
echo "  M:       ${M}"
echo "  delta:   ${DELTA}"
echo "  Condition tag: ${COND_TAG}"
echo "============================================================"

mkdir -p "${RESULTS_DIR}/data" \
         "${RESULTS_DIR}/CR"   \
         "${RESULTS_DIR}/GLM"  \
         "${RESULTS_DIR}/Q2"   \
         "${RESULTS_DIR}/LM"   \
         "${LOG_DIR}"

# ============================================================
# STEP 1 — GENERATE SIMULATED DATASETS
# Array: task index = seed number (1 .. N_SEEDS)
# ============================================================
GEN_JOB=$(sbatch \
  --parsable \
  --job-name="qsr_gen_${COND_TAG}" \
  --array="1-${N_SEEDS}" \
  -q "${QUEUE}" \
  --mem="${MEM_DATA}" \
  --time="${TIME_DATA}" \
  -N 1 \
  -n "${CPUS}" \
  --output="${LOG_DIR}/gen_${COND_TAG}_%a.out" \
  --error="${LOG_DIR}/gen_${COND_TAG}_%a.err" \
  --mail-type=FAIL \
  --mail-user="${MAIL_USER}" \
  --wrap="
module unload gnu7; module load gnu9 R

echo \"========================================\"
echo \"Job ID:      \${SLURM_JOB_ID}\"
echo \"Array Task:  \${SLURM_ARRAY_TASK_ID}\"
echo \"Node:        \$(hostname)\"
echo \"Node List:   \${SLURM_JOB_NODELIST}\"
echo \"CPU Model:   \$(lscpu | grep 'Model name' | sed 's/Model name:[[:space:]]*//')\"
echo \"Start Time:  \$(date)\"
echo \"========================================\"

echo \"Generating seed \${SLURM_ARRAY_TASK_ID} at: \$(date)\"
${RSCRIPT} ${SCRIPT_DIR}/simulate_data_quasar2.R \
  \${SLURM_ARRAY_TASK_ID} \
  ${RESULTS_DIR}/data \
  ${N_LO} ${N_HI} ${M} ${DELTA}
echo \"Done at: \$(date)\"
")

echo "Step 1 submitted: job ${GEN_JOB} (array 1-${N_SEEDS})"

# ============================================================
# STEP 2 — FIT EACH METHOD  (depends on Step 1)
# One array per method; each task picks its seed's RDS file by
# sorting the data dir and selecting the Nth entry.
# ============================================================

DATA_DIR="${RESULTS_DIR}/data"

submit_method_array() {
  local method_tag="$1"
  local r_script="$2"
  local out_subdir="${RESULTS_DIR}/${method_tag}"

  sbatch \
    --parsable \
    --job-name="qsr_${method_tag}_${COND_TAG}" \
    --array="1-${N_SEEDS}" \
    --dependency="afterok:${GEN_JOB}" \
    -q "${QUEUE}" \
    --mem="${MEM_METHOD}" \
    --time="${TIME_METHOD}" \
    -N 1 \
    -n "${CPUS}" \
    --output="${LOG_DIR}/${method_tag}_${COND_TAG}_%a.out" \
    --error="${LOG_DIR}/${method_tag}_${COND_TAG}_%a.err" \
    --mail-type=FAIL \
    --mail-user="${MAIL_USER}" \
    --wrap="
module unload gnu7; module load gnu9 R

echo \"========================================\"
echo \"Job ID:      \${SLURM_JOB_ID}\"
echo \"Array Task:  \${SLURM_ARRAY_TASK_ID}\"
echo \"Method:      ${method_tag}\"
echo \"Node:        \$(hostname)\"
echo \"Node List:   \${SLURM_JOB_NODELIST}\"
echo \"CPU Model:   \$(lscpu | grep 'Model name' | sed 's/Model name:[[:space:]]*//')\"
echo \"Start Time:  \$(date)\"
echo \"========================================\"

RDS_FILE=${RESULTS_DIR}/data/sim_data_seed\$(printf '%02d' \${SLURM_ARRAY_TASK_ID})_N${N_LO}-${N_HI}_M${M}_d${DELTA}.rds
if [[ ! -f \"\${RDS_FILE}\" ]]; then
  echo \"ERROR: expected file not found: \${RDS_FILE}\"
  exit 1
fi
echo \"[${method_tag}] Task \${SLURM_ARRAY_TASK_ID}: \${RDS_FILE} at: \$(date)\"
${RSCRIPT} ${r_script} \"\${RDS_FILE}\" \"${out_subdir}\"
echo \"[${method_tag}] Done at: \$(date)\"
"
}

CR_JOB=$(submit_method_array  "CR"  "${SCRIPT_DIR}/run_fit_QuASAR2CR.R")
GLM_JOB=$(submit_method_array "GLM" "${SCRIPT_DIR}/run_fit_QuASAR2_GLM.R")
Q2_JOB=$(submit_method_array  "Q2"  "${SCRIPT_DIR}/run_fit_QuASAR2.R")
LM_JOB=$(submit_method_array  "LM"  "${SCRIPT_DIR}/run_fit_LM.R")

echo "Step 2 submitted:"
echo "  QuASAR2CR  : job ${CR_JOB}"
echo "  QuASAR_GLM : job ${GLM_JOB}"
echo "  QuASAR2    : job ${Q2_JOB}"
echo "  LM         : job ${LM_JOB}"

# ============================================================
# STEP 3 — AGGREGATE + PLOT  (depends on all Step 2 arrays)
# ============================================================
OUT_PDF="${RESULTS_DIR}/QuASAR2_power_analysis_${COND_TAG}.pdf"

AGG_JOB=$(sbatch \
  --parsable \
  --job-name="qsr_agg_${COND_TAG}" \
  --dependency="afterok:${CR_JOB}:${GLM_JOB}:${Q2_JOB}:${LM_JOB}" \
  -q "${QUEUE}" \
  --mem="${MEM_AGG}" \
  --time="${TIME_AGG}" \
  -N 1 \
  -n "${CPUS}" \
  --output="${LOG_DIR}/agg_${COND_TAG}.out" \
  --error="${LOG_DIR}/agg_${COND_TAG}.err" \
  --mail-type=END,FAIL \
  --mail-user="${MAIL_USER}" \
  --wrap="
module unload gnu7; module load gnu9 R

echo \"========================================\"
echo \"Job ID:      \${SLURM_JOB_ID}\"
echo \"Node:        \$(hostname)\"
echo \"Node List:   \${SLURM_JOB_NODELIST}\"
echo \"CPU Model:   \$(lscpu | grep 'Model name' | sed 's/Model name:[[:space:]]*//')\"
echo \"Start Time:  \$(date)\"
echo \"========================================\"

echo \"Aggregating results at: \$(date)\"
${RSCRIPT} ${SCRIPT_DIR}/simulation_result_analyses.R \
  ${RESULTS_DIR} \
  ${OUT_PDF}
echo \"Finished at: \$(date)\"
")

echo "Step 3 submitted: job ${AGG_JOB}"
echo ""
echo "All jobs queued. Monitor with: squeue -u \$USER"
echo "Final plot will be: ${OUT_PDF}"
echo "Finished submitting at: $(date)"

# ============================================================
# EXAMPLES — sweeping conditions
# ============================================================
# To sweep across multiple conditions, submit this script in a loop:
#
#   for M in 20 100 500; do
#     for DELTA in 0.05 0.10 0.20; do
#       sbatch run_simulations.sh --M $M --delta $DELTA --N_lo 60 --N_hi 300
#     done
#   done
#
# Each sbatch call produces its own job chain and its own PDF.
# To produce a combined plot across all conditions after everything finishes:
#
#   Rscript R/simulation_result_analyses.R results/ combined_power_analysis.pdf