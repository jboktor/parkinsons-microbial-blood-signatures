#!/bin/bash
#SBATCH --job-name=setup-deeprc-v2
#SBATCH --output=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/setup_deeprc_env_v2_%j.out
#SBATCH --error=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/setup_deeprc_env_v2_%j.err
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

## Remove old env
echo "=== Removing old env ==="
mamba remove -n immuneml_deeprc --all -y 2>&1 || true

## Step 1: Create env with pytorch-cuda from conda (all deps resolved together)
echo "=== Creating env with GPU PyTorch + immuneml ==="
mamba create -n immuneml_deeprc \
    -c pytorch -c nvidia -c conda-forge -c bioconda \
    --channel-priority flexible \
    python=3.11 \
    pytorch pytorch-cuda=12.1 \
    immuneml \
    -y

## Step 2: Install DeepRC (pip, no-deps to avoid breaking conda packages)
echo "=== Installing DeepRC ==="
mamba activate immuneml_deeprc
pip install --no-dependencies git+https://github.com/widmi/widis-lstm-tools.git
pip install --no-dependencies git+https://github.com/ml-jku/DeepRC.git

## Step 3: Verify
echo "=== Verifying ==="
python -c "
import torch
print('PyTorch:', torch.__version__)
print('CUDA available:', torch.cuda.is_available())
if torch.cuda.is_available():
    print('GPU:', torch.cuda.get_device_name(0))
    print('CUDA version:', torch.version.cuda)
else:
    print('WARNING: CUDA not available')
from immuneML.ml_methods.classifiers.DeepRC import DeepRC
from immuneML.encodings.deeprc.DeepRCEncoder import DeepRCEncoder
print('DeepRC imports OK')
import scipy
print('scipy OK:', scipy.__version__)
"

## Step 4: Run minimal test
echo "=== Running minimal DeepRC test ==="
RESULTS_DIR=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/minimal_test_deeprc
rm -rf ${RESULTS_DIR}/*
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

immune-ml /resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/ZZZ_minimal_test_deeprc.yaml .

echo "=== Done ==="
