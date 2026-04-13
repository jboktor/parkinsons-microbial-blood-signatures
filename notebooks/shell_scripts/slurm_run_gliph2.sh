#!/bin/bash
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=15G
#SBATCH --time=7-00:00:00   # time limit
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL   # Notify on failure.
#SBATCH --output=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/gliph2/gliph2_tcrb_%j.log
#SBATCH --error=/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/gliph2/gliph2_tcrb_%j.err

source /home/${USER}/.bashrc

cd /resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/gliph2/tcrb_results_2025-07-16

/resnick/groups/MazmanianLab/jboktor/software/gliph2/irtools.centos -c /resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/input/gliph2/param_files/gliph2_tcrb.cfg
