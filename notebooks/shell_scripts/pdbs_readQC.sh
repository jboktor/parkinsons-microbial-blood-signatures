#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=4:00:00   # walltime
#SBATCH --ntasks=1 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=20G   # memory per CPU core
#SBATCH -J "ReadQC"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=<jboktor>@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output=/central/scratch/jbok/slurmdump/ReadQC_%j.out


while getopts r:s:o: option
do
case "${option}"
in
r) ROOTDIR=${OPTARG};;
s) SAMPLEID=${OPTARG};;
o) OUTPUTDIR=${OPTARG};;
esac
done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

source /home/${USER}/.bashrc
source activate pdbm
cd $OUTPUTDIR

echo "PROCESSING SAMPLE: "$SAMPLEID
SAMPLE_IN="${ROOTDIR}$SAMPLEID"
SAMPLE_OUT="${OUTPUTDIR}$SAMPLEID"

# This command conducts the following:
# ______________________________________
# trims low quality reads from both ends with PHRED  Q<10
# # qtrim=rl trimq=10
# after quality trimming, this will calcluate the average quality and 
# remove reads with an average PHRED Q<10
# # maq=10
# This will filter out reads that have an average entropy of under 0.6. 
# with a sliding window of 50 and k-mer length of five
# This is low read complexity filter
# # entropy=0.5 entropywindow=50 entropyk=5

# STDOUT = "bbduk_stdout_$SAMPLEOUT.txt"

bbduk.sh qtrim=rl trimq=10 maq=10 entropy=0.5 entropywindow=50 entropyk=5 overwrite=true \
    in1=${SAMPLE_IN}_1.fq.gz \
    in2=${SAMPLE_IN}_2.fq.gz \
    out1=${SAMPLE_OUT}_1.fq.gz \
    out2=${SAMPLE_OUT}_2.fq.gz >& "${SAMPLEID}_bbduk_stdout.txt"
    