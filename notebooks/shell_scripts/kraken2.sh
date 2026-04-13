#!/bin/bash 

# This script takes in only the name of sample
# not including forward/reverse/singleton ext ie (Sample, not Sample_2.fq.gz)
# and concatenates F/R/S reads and uses as input for kraken2 classification using 
# three different RefDBs (RefSeqPlusPF, WoL, and UHGG)

while getopts s: option
do
case "${option}"
in
s) SAMPLEID=${OPTARG};;
esac
done

#______________________________________________________________________________
#                     Slurm Construction Section
#______________________________________________________________________________

#SBATCH --time=10:00:00   # walltime
#SBATCH --ntasks=4 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=75G   # memory per CPU core
#SBATCH -J "${SAMPLEID}_Krak2"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output=/resnick/scratch/jbok/kraken2_stdout/${SAMPLEID}_Kraken2_%j.out

#______________________________________________________________________________

source /home/${USER}/.bashrc
INPUTDIR='/resnick/groups/MazmanianLab/jboktor/PDMBS/workflow/WGS/clean_fastqs/'
OUTPUTDIR='/resnick/groups/MazmanianLab/jboktor/PDMBS/workflow/WGS/results/kraken2/' 
# SEQOUTDIR='/resnick/scratch/jbok/kraken_sequences/' 

SAMPLEINPUT="${INPUTDIR}${SAMPLEID}"
TEMPDIR='/resnick/scratch/jbok/krakenScratch/'
READSDIR="${TEMPDIR}${SAMPLEID}.fq.gz"

echo "PROCESSING SAMPLE: "$SAMPLEID
echo "SAMPLE INPUT: "$SAMPLEINPUT
echo "SAMPLE OUTPUT: "$READSDIR

# make a concatenated fastq file with forward/reverse/singleton reads
cat ${SAMPLEINPUT}"_1.fq.gz" ${SAMPLEINPUT}"_2.fq.gz" ${SAMPLEINPUT}"_single.fq.gz" > $READSDIR
#______________________________________________________________________________

refseq_kraken_run="kraken2 --db /resnick/groups/MazmanianLab/jboktor/Downloads/refseq_pluspf_v4/ \
--threads 4 \
--gzip-compressed \
--classified-out "$OUTPUTDIR'RefSeqPlusPF_mapped/'${SAMPLEID}'__classified_RefSeqPlusPF.fastq'" \
--report "$OUTPUTDIR'RefSeqPlusPF_mapped/'${SAMPLEID}'__report_RefSeqPlusPF.tsv'" ${READSDIR}"
echo $refseq_kraken_run
$refseq_kraken_run
gzip $OUTPUTDIR"RefSeqPlusPF_mapped/"${SAMPLEID}"__classified_RefSeqPlusPF.fastq"


uhgg_kraken_run="kraken2 --db /resnick/groups/MazmanianLab/jboktor/Downloads/uhgg_kraken2-db/ \
--threads 4 \
--gzip-compressed \
--classified-out "$OUTPUTDIR'UHGG_mapped/'${SAMPLEID}'__classified_UHGG.fastq'" \
--report "$OUTPUTDIR'UHGG_mapped/'${SAMPLEID}'__report_UHGG.tsv'"  ${READSDIR}"
echo $uhgg_kraken_run
$uhgg_kraken_run
gzip $OUTPUTDIR"UHGG_mapped/"${SAMPLEID}"__classified_UHGG.fastq"


WoL_kraken_run="kraken2 --db /resnick/groups/MazmanianLab/jboktor/WebOfLife/databases/kraken2/ \
--threads 4 \
--gzip-compressed \
--classified-out "$OUTPUTDIR'WoL_mapped/'${SAMPLEID}'__classified_WoL.fastq'" \
--report "$OUTPUTDIR'WoL_mapped/'${SAMPLEID}'__report_WoL.tsv'" ${READSDIR}"
echo $WoL_kraken_run
$WoL_kraken_run
gzip $OUTPUTDIR"WoL_mapped/"${SAMPLEID}"__classified_WoL.fastq"

rm $READSDIR
