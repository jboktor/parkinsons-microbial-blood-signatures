#!/bin/bash
#SBATCH --job-name=deeprc_v2
#SBATCH --output=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/deeprc_v2_%A_%a.out
#SBATCH --error=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/deeprc_v2_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_h200:1

source /home/${USER}/.bashrc
module load cuda/12.2.1-gcc-11.3.1-sdqrj2e

MANIFEST=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/deeprc_v2_train_manifest.txt
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $MANIFEST)
CONFIG_NAME=$(echo "$LINE" | awk '{print $1}')
SEED=$(echo "$LINE" | awk '{print $2}')

CONFIGS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/deeprc_v2
RESULTS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/deeprc_v2

RESULTS_DIR=${RESULTS}/${CONFIG_NAME}_seed${SEED}

echo "=== ${CONFIG_NAME} seed=${SEED} ==="
echo "Start: $(date)"
nvidia-smi

mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

mamba activate immuneml_deeprc

python -c "import torch; print('CUDA:', torch.cuda.is_available(), '| GPU:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none')"

# Inject the seed into the YAML by creating a seeded copy
TMP_YAML=/tmp/deeprc_v2_${CONFIG_NAME}_seed${SEED}.yaml
cp ${CONFIGS}/${CONFIG_NAME}.yaml ${TMP_YAML}

# ImmuneML doesn't have an explicit seed arg — use PYTHONHASHSEED and torch seed via env
export PYTHONHASHSEED=${SEED}

immune-ml ${TMP_YAML} .

echo "Exit: $?"
echo "End: $(date)"
