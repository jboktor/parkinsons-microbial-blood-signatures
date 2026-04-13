#!/bin/bash
#SBATCH --job-name=minimal-deeprc-test
#SBATCH --output=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/minimal_deeprc_%j.out
#SBATCH --error=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/minimal_deeprc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=END,FAIL

## Initialize work environment
source /home/${USER}/.bashrc
module load cuda/12.2.1-gcc-11.3.1-sdqrj2e

export TORCH_HOME=/resnick/groups/MazmanianLab/jboktor/cache/torch
export TRANSFORMERS_CACHE=/resnick/groups/MazmanianLab/jboktor/cache

## Use scratch for tmp
export TMPDIR=/resnick/scratch/jbok/deeprc_test/${SLURM_JOBID}
mkdir -p ${TMPDIR}

## Setup results directory
RESULTS_DIR=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/minimal_test_deeprc
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

## Run ImmuneML DeepRC (separate env with conda GPU PyTorch + deeprc)
mamba activate immuneml_deeprc
echo "PyTorch CUDA check:"
python -c "import torch; print('CUDA:', torch.cuda.is_available(), '| device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none')"

echo "=== GPU info ==="
nvidia-smi
echo "=== Starting ImmuneML DeepRC ==="
echo "Start time: $(date)"

immune-ml /resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/ZZZ_minimal_test_deeprc.yaml .

echo "End time: $(date)"
echo "=== Done ==="
