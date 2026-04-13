#!/bin/bash
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=24
#SBATCH --mem=1000GB
#SBATCH --time=7-00:00:00   # time limit
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL   # Notify on failure.
#SBATCH --output=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/giana/giana_tcrb_%j.log
#SBATCH --error=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/giana/giana_tcrb_%j.err

source /home/${USER}/.bashrc

source activate giana

python /resnick/groups/MazmanianLab/jboktor/git/GIANA/GIANA4.1.py \
    -f /resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/giana/input/full_tcrb_2025-07-10.tsv \
    -S 3 \
    -t 10 \
    -N 24 \
    -v \
    -o /resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/giana/results
