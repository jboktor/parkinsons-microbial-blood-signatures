#!/bin/bash
#SBATCH --job-name=S2S_deeprc
#SBATCH --output=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/S2S_test_deeprc_%j.out
#SBATCH --error=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/S2S_test_deeprc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:h100:1

## Initialize
source /home/${USER}/.bashrc
module load cuda/12.2.1-gcc-11.3.1-sdqrj2e

CONFIGS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs
RESULTS=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results

echo "=== Running S2S_test_deeprc (GPU) ==="
echo "Start time: $(date)"
nvidia-smi

RESULTS_DIR=${RESULTS}/S2S_test_deeprc
rm -rf ${RESULTS_DIR}
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

mamba activate immuneml_deeprc

echo "PyTorch CUDA check:"
python -c "import torch; print('CUDA:', torch.cuda.is_available(), '| device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none')"

# Use cuda:0 version of the config
sed 's/pytorch_device_name: cpu.*$/pytorch_device_name: cuda:0/' ${CONFIGS}/S2S_test_deeprc.yaml > /tmp/S2S_test_deeprc_gpu.yaml

immune-ml /tmp/S2S_test_deeprc_gpu.yaml .

echo "Exit code: $?"
echo "End time: $(date)"
echo "=== Done ==="
