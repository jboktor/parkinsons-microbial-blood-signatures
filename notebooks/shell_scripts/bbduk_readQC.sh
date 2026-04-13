#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=4:00:00   # walltime
#SBATCH --ntasks=1 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=20G   # memory per CPU core
#SBATCH -J "ReadQC"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output=/resnick/scratch/jbok/slurmdump/ReadQC_%j.out
#SBATCH --error=/resnick/scratch/jbok/slurmdump/ReadQC_%j.err

# This scripts concuts bbduk read quality trimming on a per-file basis
# forward, reverse, and singleton reads are processed independently
# This script requires an input directory, a gziped fastq file, 
# and an output directory 
# ie.
# sbatch bbduk_readQC.sh \
# -i path/to/rawfastqfiles \
# -s SY-PDYV487GW8_single.fq.gz |
# -o path/to/cleanfastqfiles


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
cd /resnick/groups/MazmanianLab/jboktor/PDMBS/workflow/WGS/clean_fastqs_stats/

# This command conducts the following:
# ______________________________________
# trims low quality reads from both ends with PHRED  Q<10
# # qtrim=rl trimq=10
# after quality trimming, this will calcluate the average quality and 
# remove reads with an average PHRED Q<10
# # maq=10

bbduk.sh qtrim=rl trimq=10 maq=10 overwrite=true -Xmx15g \
    in=${SAMPLE_IN} \
    out=${SAMPLE_OUT} >& "${SAMPLE_NAME}_bbduk_stdout.txt"
