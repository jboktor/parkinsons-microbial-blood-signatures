#!/bin/bash
#SBATCH --job-name=fix-deeprc-cuda
#SBATCH --output=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/fix_deeprc_cuda_%j.out
#SBATCH --error=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/fix_deeprc_cuda_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

source /home/${USER}/.bashrc
module load cuda/12.2.1-gcc-11.3.1-sdqrj2e

echo "=== GPU info ==="
nvidia-smi

echo "=== Step 1: Check current PyTorch ==="
mamba activate immuneml_deeprc
python -c "import torch; print('Before:', torch.__version__, '| CUDA:', torch.cuda.is_available())"

echo "=== Step 2: Install CUDA PyTorch via pip (keeping all other conda packages) ==="
# Use --index-url (not --extra-index-url) to ensure we get the cu121 build
# Only replace torch, torchvision, torchaudio — leave everything else alone
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

echo "=== Step 3: Verify CUDA PyTorch ==="
python -c "
import torch
print('After:', torch.__version__)
print('CUDA available:', torch.cuda.is_available())
if torch.cuda.is_available():
    print('GPU:', torch.cuda.get_device_name(0))
    print('CUDA version:', torch.version.cuda)
else:
    print('CUDA NOT available - check error above')
"

echo "=== Step 4: Verify immuneML + DeepRC still import ==="
python -c "
from immuneML.ml_methods.classifiers.DeepRC import DeepRC
from immuneML.encodings.deeprc.DeepRCEncoder import DeepRCEncoder
import scipy
print('DeepRC imports OK')
print('scipy:', scipy.__version__)
"

echo "=== Step 5: Run minimal DeepRC test with cuda:0 ==="
# Temporarily patch the YAML to use cuda:0
YAML=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/ZZZ_minimal_test_deeprc.yaml
RESULTS_DIR=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/minimal_test_deeprc_gpu
rm -rf ${RESULTS_DIR}
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

# Create a temp YAML with cuda:0
sed 's/pytorch_device_name: cpu.*$/pytorch_device_name: cuda:0/' ${YAML} > /tmp/deeprc_gpu_test.yaml

immune-ml /tmp/deeprc_gpu_test.yaml .

echo "=== Done ==="
