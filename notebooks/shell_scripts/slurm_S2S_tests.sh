#!/bin/bash
# Master script to submit all 6 S2S transfer-study ML tests to SLURM
# Usage: bash slurm_S2S_tests.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIGS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs
RESULTS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results
LOG_DIR=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs

# CPU models (immuneml env)
for model in kmer_logreg kmer_svm kmer_rf kmer_abundance pubclone; do
    echo "Submitting S2S_test_${model}..."
    sbatch ${SCRIPT_DIR}/slurm_S2S_cpu.sh ${model}
done

# GPU model (immuneml_deeprc env)
echo "Submitting S2S_test_deeprc..."
sbatch ${SCRIPT_DIR}/slurm_S2S_deeprc.sh
