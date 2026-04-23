#!/bin/bash
#SBATCH --job-name=compairr_rerun
#SBATCH --output=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/compairr_rerun_%A_%a.out
#SBATCH --error=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/compairr_rerun_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=900GB
#SBATCH --time=24:00:00
#SBATCH --partition=expansion

source /home/${USER}/.bashrc

MANIFEST=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/compairr_rerun_manifest.txt
CONFIG_NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $MANIFEST)

CONFIGS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/compairr_rerun
RESULTS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/compairr_rerun

echo "=== ${CONFIG_NAME} (CompAIRR train+apply, 900GB, ignore_genes=true) ==="
echo "Start: $(date)"

RESULTS_DIR=${RESULTS}/${CONFIG_NAME}
rm -rf ${RESULTS_DIR}
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

mamba activate immuneml
immune-ml ${CONFIGS}/${CONFIG_NAME}.yaml .

echo "Exit: $?"
echo "End: $(date)"
