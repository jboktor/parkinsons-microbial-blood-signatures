#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --partition=expansion

## Initialize
source /home/${USER}/.bashrc

# Args: $1 = design (S2S or LOSO), $2 = yaml filename (without .yaml)
DESIGN=$1
CONFIG_NAME=$2
CONFIGS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/${DESIGN}
RESULTS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/${DESIGN}

echo "=== Running ${DESIGN}/${CONFIG_NAME} ==="
echo "Start time: $(date)"

RESULTS_DIR=${RESULTS}/${CONFIG_NAME}
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

mamba activate immuneml
immune-ml ${CONFIGS}/${CONFIG_NAME}.yaml .

echo "Exit code: $?"
echo "End time: $(date)"
echo "=== Done ==="
