#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=1-00:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=100G   # memory per CPU core
#SBATCH -J "Kraken2-testrun"   # job name
#SBATCH --mail-user=<jboktor>@caltech.edu   # email address

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     User Construction Section

source /home/${USER}/.bashrc
source activate wol

kraken2 --db /resnick/groups/MazmanianLab/jboktor/WebOfLife/databases/kraken2/ \
-threads 2 \
--gzip-compressed \
--classified-out /resnick/groups/MazmanianLab/jboktor/PDBM/WoL_mapped/BF-1002__classified_WoL.tsv \
--unclassified-out /resnick/groups/MazmanianLab/jboktor/PDBM/WoL_mapped/BF-1002__unclassified_WoL.tsv \
--report /resnick/groups/MazmanianLab/jboktor/PDBM/WoL_mapped/BF-1002__report_WoL.tsv \
<(cat /resnick/groups/MazmanianLab/jboktor/PDBM/test_input/BF-1002_single.fq.gz && /home/jboktor/bbmap/bbmerge-auto.sh in1=/resnick/groups/MazmanianLab/jboktor/PDBM/test_input/BF-1002_1.fq.gz in2=/resnick/groups/MazmanianLab/jboktor/PDBM/test_input/BF-1002_2.fq.gz)

# kraken2 --db /resnick/groups/MazmanianLab/jboktor/Downloads/refseq_pluspf_v4/ \
# --threads 4 \
# --gzip-compressed \
# --classified-out testrun/test_classfied_refseq.tsv \
# --unclassified-out testrun/test_unclassified_refseq.tsv \
# --report testrun/test_report_refseq.tsv \
# <(cat testdatacp/BF-1003_single.fq.gz && $HOME/bbmap/bbmerge-auto.sh in1=testdatacp/BF-1003_1.fq.gz in2=testdatacp/BF-1003_2.fq.gz) 


# <(cat testdatacp/BF-1003_single.fq.gz && $HOME/FLASH-1.2.11-Linux-x86_64/flash --max-overlap 150  testdatacp/BF-1003_1.fq.gz testdatacp/BF-1003_2.fq.gz) 


# kraken2 --db /resnick/groups/MazmanianLab/jboktor/Downloads/uhgg_kraken2-db/ \
# --threads 16 \
# --gzip-compressed \
# --classified-out test_classified_uhgg.tsv \
# --unclassified-out test_unclassified_uhgg.tsv \
# --report test_report_uhgg \
# testdata/618fa114-17bc-49fc-95c1-56667e68cebc_cram_to_fastq_workflow_8d71236f-cdd5-4d19-9c62-75f5484601d1_call-cram_to_fastq_attempt-2_BF-1003_single.fq.gz

#--db /resnick/groups/MazmanianLab/jboktor/WebOfLife/databases/kraken2/ \
#--db /resnick/groups/MazmanianLab/jboktor/Downloads/refseq_pluspf_v4/ \