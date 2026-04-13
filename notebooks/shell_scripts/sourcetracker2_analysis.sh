#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=6:00:00   # walltime
#SBATCH --ntasks=10 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32G   # memory per CPU core
#SBATCH -J "Sourcetracker2"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=<jboktor>@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output=/resnick/scratch/jbok/slurmdump/Sourcetracker2_%j.out


source /home/${USER}/.bashrc
source activate st2
cd /resnick/groups/MazmanianLab/jboktor/PDBM


sourcetracker2 gibbs -i classification/RefSeqPlusPF_kraken2.biom -m RefSeqPlusPF_mapping.txt -o st2analysis/ --jobs 10 
