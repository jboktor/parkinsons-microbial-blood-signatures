#!/bin/bash
# Launch all S2S and LOSO ImmuneML jobs on SLURM
# Usage: bash launch_full_study.sh [S2S|LOSO|all]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIGS=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs
RESULTS=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results
LOG_DIR=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs

MODE=${1:-all}
SUBMITTED=0

submit_job() {
    local design=$1
    local config_name=$2
    sbatch \
        --job-name="${design}_${config_name}" \
        --output="${LOG_DIR}/${design}_${config_name}_%j.out" \
        --error="${LOG_DIR}/${design}_${config_name}_%j.err" \
        ${SCRIPT_DIR}/immuneml_slurm_full_study_cpu.sh ${design} ${config_name}
    SUBMITTED=$((SUBMITTED + 1))
}

if [ "$MODE" = "S2S" ] || [ "$MODE" = "all" ]; then
    echo "=== Submitting S2S jobs ==="
    for yaml in ${CONFIGS}/S2S/S2S_*.yaml; do
        config_name=$(basename $yaml .yaml)
        submit_job "S2S" "$config_name"
    done
fi

if [ "$MODE" = "LOSO" ] || [ "$MODE" = "all" ]; then
    echo "=== Submitting LOSO jobs ==="
    for yaml in ${CONFIGS}/LOSO/LOSO_*.yaml; do
        config_name=$(basename $yaml .yaml)
        submit_job "LOSO" "$config_name"
    done
fi

echo "=== Submitted ${SUBMITTED} jobs ==="
