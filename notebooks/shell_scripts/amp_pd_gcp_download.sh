#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=1 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1000G   # memory per CPU core
#SBATCH -J "GCP_Download"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=<jboktor>@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output=/resnick/scratch/jbok/slurmdump/compile_reads_%j.out


while getopts f: option
do
case "${option}"
in
f) GCP_FOLDER=${OPTARG};;
esac
done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cd /resnick/groups/MazmanianLab/jboktor/PDBM/workflow_downloads/WGS

gsutil -u amp-pd-gcp-joeb -m cp -n -r ${GCP_FOLDER} . 

