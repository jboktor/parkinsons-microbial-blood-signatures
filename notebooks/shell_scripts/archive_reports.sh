#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=1 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=100G   # memory per CPU core
#SBATCH -J "Archive_Krakens"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=<jboktor>@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output=/central/scratch/jbok/slurmdump/Archive_Krakens_%j.out

source /home/${USER}/.bashrc
source activate wol

cd /central/groups/MazmanianLab/joeB/PDMBS/classification
today=$(date +"%Y-%m-%d")

# place archived files in scratch folder for now (remove scratch prefix later)
tar -cvf /central/scratch/jbok/RefSeqPlusPF_mapped_reports_${today}.tar RefSeqPlusPF_mapped/*report_RefSeqPlusPF.tsv
tar -cvf /central/scratch/jbok/UHGG_mapped_reports_${today}.tar UHGG_mapped/*report_UHGG.tsv
tar -cvf /central/scratch/jbok/WoL_mapped_reports_${today}.tar WoL_mapped/*report_WoL.tsv


# kraken-biom RefSeqPlusPF_mapped/*_report_RefSeqPlusPF.tsv \
# -o RefSeqPlusPF_kraken2_${today}.biom \
# --gzip --max D --min S 

# kraken-biom UHGG_mapped/*_report_UHGG.tsv \
# -o UHGG_kraken2_${today}.biom \
# --gzip --max D --min S 

# kraken-biom WoL_mapped/*_report_WoL.tsv \
# -o WoL_kraken2_${today}.biom \
# --gzip --max D --min S 

