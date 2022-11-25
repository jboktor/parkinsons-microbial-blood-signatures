#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=4:00:00   # walltime
#SBATCH --ntasks=1 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=50G   # memory per CPU core
#SBATCH -J "ReadQC"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=<jboktor>@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output=/central/scratch/jbok/slurmdump/ReadQC_%j.out


while getopts i:s:o: option
do
case "${option}"
in
i) INPUTDIR=${OPTARG};;
s) SAMPLEID=${OPTARG};;
o) OUTPUTDIR=${OPTARG};;
esac
done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

source /home/${USER}/.bashrc
source activate pdbm

SAMPLE_IN="${INPUTDIR}$SAMPLEID"
SAMPLE_OUT="${OUTPUTDIR}$SAMPLEID"
SAMPLE_NAME=`echo ${SAMPLEID} | sed 's/.fq.gz//'`

echo "PROCESSING SAMPLE: "$SAMPLEID

# move to stats output directory 
cd /central/groups/MazmanianLab/joeB/PDMBS/workflow/WGS/clean_fastqs_stats/

# This command conducts the following:
# ______________________________________
# trims low quality reads from both ends with PHRED  Q<10
# # qtrim=rl trimq=10
# after quality trimming, this will calcluate the average quality and 
# remove reads with an average PHRED Q<10
# # maq=10

bbduk.sh qtrim=rl trimq=10 maq=10 overwrite=true \
    in=${SAMPLE_IN} \
    out=${SAMPLE_OUT} >& "${SAMPLE_NAME}_bbduk_stdout.txt"

# bbduk.sh qtrim=rl trimq=10 maq=10 overwrite=true \
#     in1=${SAMPLE_IN}_1.fq.gz \
#     in2=${SAMPLE_IN}_2.fq.gz \
#     out1=${SAMPLE_OUT}_1.fq.gz \
#     out2=${SAMPLE_OUT}_2.fq.gz >& "${SAMPLEID}_bbduk_stdout.txt"