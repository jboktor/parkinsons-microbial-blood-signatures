#!/bin/bash
#SBATCH --job-name=minimal-deeprc-test
#SBATCH --output=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/minimal_deeprc_%j.out
#SBATCH --error=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/.cluster_runs/minimal_deeprc_%j.err
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

export TORCH_HOME=/central/groups/MazmanianLab/jboktor/cache/torch
export TRANSFORMERS_CACHE=/central/groups/MazmanianLab/jboktor/cache

## Use scratch for tmp
export TMPDIR=/central/scratch/jbok/deeprc_test/${SLURM_JOBID}
mkdir -p ${TMPDIR}

## Setup results directory
RESULTS_DIR=/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/results/minimal_test_deeprc
mkdir -p ${RESULTS_DIR}
cd ${RESULTS_DIR}

## Run ImmuneML DeepRC
mamba activate immuneml

echo "=== GPU info ==="
nvidia-smi
echo "=== Starting ImmuneML DeepRC ==="
echo "Start time: $(date)"

immune-ml /central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs/ZZZ_minimal_test_deeprc.yaml .

echo "End time: $(date)"
echo "=== Done ==="
