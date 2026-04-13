#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:nvidia_h200:1

## Initialize
source /home/${USER}/.bashrc
module load cuda/12.2.1-gcc-11.3.1-sdqrj2e

# Args: $1 = design (S2S or LOSO), $2 = yaml filename (without .yaml)
DESIGN=$1
CONFIG_NAME=$2
CONFIGS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/${DESIGN}
RESULTS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/${DESIGN}

echo "=== Running ${DESIGN}/${CONFIG_NAME} (GPU DeepRC) ==="
echo "Start time: $(date)"
nvidia-smi

RESULTS_DIR=${RESULTS}/${CONFIG_NAME}
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

mamba activate immuneml_deeprc

echo "PyTorch CUDA check:"
python -c "import torch; print('CUDA:', torch.cuda.is_available(), '| device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none')"

immune-ml ${CONFIGS}/${CONFIG_NAME}.yaml .

echo "Exit code: $?"
echo "End time: $(date)"
echo "=== Done ==="
