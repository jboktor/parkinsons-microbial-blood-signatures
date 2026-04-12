#!/bin/bash
#SBATCH --job-name=setup-deeprc-env
#SBATCH --output=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/setup_deeprc_env_%j.out
#SBATCH --error=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/setup_deeprc_env_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

## Initialize work environment
source /home/${USER}/.bashrc
module load cuda/12.2.1-gcc-11.3.1-sdqrj2e

echo "=== GPU info ==="
nvidia-smi
echo "=== CUDA version ==="
nvcc --version 2>/dev/null || echo "nvcc not found, using module cuda"

## Remove old env if it exists
echo "=== Removing old immuneml_deeprc env ==="
mamba remove -n immuneml_deeprc --all -y 2>&1 || true

## Create fresh env with conda immuneml
echo "=== Creating immuneml_deeprc env ==="
mamba create -n immuneml_deeprc -c conda-forge -c bioconda immuneml python=3.11 -y

## Activate and install GPU packages
echo "=== Installing GPU PyTorch ==="
mamba activate immuneml_deeprc
pip install --force-reinstall torch --extra-index-url https://download.pytorch.org/whl/cu121

echo "=== Installing DeepRC + deps ==="
pip install --no-dependencies git+https://github.com/widmi/widis-lstm-tools.git
pip install --no-dependencies git+https://github.com/ml-jku/DeepRC.git
pip install tensorboard h5py requests

echo "=== Verifying installation ==="
python -c "
import torch
print('PyTorch version:', torch.__version__)
print('CUDA available:', torch.cuda.is_available())
if torch.cuda.is_available():
    print('GPU:', torch.cuda.get_device_name(0))
    print('CUDA version:', torch.version.cuda)
from immuneML.ml_methods.classifiers.DeepRC import DeepRC
from immuneML.encodings.deeprc.DeepRCEncoder import DeepRCEncoder
print('DeepRC imports OK')
"

echo "=== Now running minimal DeepRC test ==="
RESULTS_DIR=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/minimal_test_deeprc
rm -rf ${RESULTS_DIR}/*
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

immune-ml /central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/ZZZ_minimal_test_deeprc.yaml .

echo "=== Done ==="
