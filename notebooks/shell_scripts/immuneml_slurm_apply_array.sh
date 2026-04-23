#!/bin/bash
#SBATCH --job-name=apply_v2
#SBATCH --output=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/apply_v2_%A_%a.out
#SBATCH --error=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/apply_v2_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=04:00:00
#SBATCH --partition=expansion

source /home/${USER}/.bashrc

MANIFEST=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/apply_v2_manifest.txt
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $MANIFEST)
CONFIG_NAME=$(echo "$LINE" | awk '{print $1}')

CONFIGS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/apply_v2
RESULTS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/apply_v2

echo "=== ${CONFIG_NAME} ==="
echo "Start: $(date)"

RESULTS_DIR=${RESULTS}/${CONFIG_NAME}
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

mamba activate immuneml
immune-ml ${CONFIGS}/${CONFIG_NAME}.yaml .

echo "Exit: $?"
echo "End: $(date)"
