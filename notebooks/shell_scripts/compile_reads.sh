#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=1 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=200G   # memory per CPU core
#SBATCH -J "MergeReads"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=<jboktor>@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output=/central/scratch/jbok/slurmdump/compile_reads_%j.out


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
source activate wol
cd $OUTPUTDIR

if [ ! -f $SAMPLEID".fq.gz" ]
then
    echo $SAMPLEID".fq.gz does not exist .... running script"
    echo "PROCESSING SAMPLE: "$SAMPLEID
    SAMPLEROOT="${ROOTDIR}$SAMPLEID"
    merge_reads="$HOME/bbmap/bbmerge-auto.sh in1=${SAMPLEROOT}"_1.fq.gz" in2=${SAMPLEROOT}"_2.fq.gz" out=${SAMPLEID}"_merged.fq.gz" outu=${SAMPLEID}"_unmerged.fq.gz""
    echo $merge_reads
    $merge_reads
    cat ${SAMPLEID}"_unmerged.fq.gz" ${SAMPLEID}"_merged.fq.gz" ${SAMPLEROOT}"_single.fq.gz" > $SAMPLEID".fq.gz"
    rm ${SAMPLEID}"_unmerged.fq.gz" ${SAMPLEID}"_merged.fq.gz" 
else
    echo "Found: "$SAMPLEID".fq.gz"
fi


# concat="cat ${SAMPLEID}"_unmerged.fq.gz" ${SAMPLEID}"_merged.fq.gz" ${SAMPLEROOT}"_single.fq.gz" > $SAMPLEID"_NEW.fq.gz""
# echo $concat
