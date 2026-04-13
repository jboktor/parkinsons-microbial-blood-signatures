#!/bin/bash

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL

while getopts "f:p:" opt
do
    case "$opt" in
        f ) FILE_NAME="$OPTARG" ;;
        p ) FILE_PATH="$OPTARG" ;;
    esac
done

# example inputs
# FILE_PATH="https://downloads.hmpdacc.org/dacc/hhs/genome/microbiome/wgs/analysis/hmwgsqc/v1/SRS024641.tar.bz2"
# FILE_NAME="SRS024641.tar.bz2"

# activate enviornment
source /home/${USER}/.bashrc
source activate pdmbsR
# move into scratch dir
cd /resnick/scratch/jbok/tmp

wget --no-check-certificate "${FILE_PATH}" && \
gsutil cp ${FILE_NAME} "gs://shotgun-metagenomic-data/HMP2/${FILE_NAME}" && \
rm "${FILE_NAME}"

