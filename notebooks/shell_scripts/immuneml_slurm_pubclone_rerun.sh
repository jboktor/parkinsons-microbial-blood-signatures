#!/bin/bash
#SBATCH --job-name=pubclone_rerun
#SBATCH --output=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/pubclone_rerun_%A_%a.out
#SBATCH --error=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/pubclone_rerun_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128GB
#SBATCH --time=08:00:00
#SBATCH --partition=expansion

source /home/${USER}/.bashrc

MANIFEST=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/pubclone_rerun_manifest.txt
CONFIG_NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $MANIFEST)

CONFIGS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/apply_v2
RESULTS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/apply_v2

echo "=== ${CONFIG_NAME} (pubclone rerun with patched env) ==="
echo "Start: $(date)"

RESULTS_DIR=${RESULTS}/${CONFIG_NAME}
rm -rf ${RESULTS_DIR}
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

mamba activate immuneml
immune-ml ${CONFIGS}/${CONFIG_NAME}.yaml .

echo "Exit: $?"
echo "End: $(date)"
