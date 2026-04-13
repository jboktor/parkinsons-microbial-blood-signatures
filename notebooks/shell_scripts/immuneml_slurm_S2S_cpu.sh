#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --partition=expansion

## Initialize
source /home/${USER}/.bashrc

MODEL_NAME=$1
CONFIGS=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs
RESULTS=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results

echo "=== Running S2S_test_${MODEL_NAME} ==="
echo "Start time: $(date)"

RESULTS_DIR=${RESULTS}/S2S_test_${MODEL_NAME}
rm -rf ${RESULTS_DIR}
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

mamba activate immuneml
immune-ml ${CONFIGS}/S2S_test_${MODEL_NAME}.yaml .

echo "Exit code: $?"
echo "End time: $(date)"
echo "=== Done ==="
