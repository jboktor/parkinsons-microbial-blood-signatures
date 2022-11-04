#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=6:00:00   # walltime
#SBATCH --ntasks=2 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=50G   # memory per CPU core
#SBATCH -J "Krak2"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=<jboktor>@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output=/central/scratch/jbok/slurmdump/%j.out

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

while getopts o:s: option
do
case "${option}"
in
o) ROOTDIR=${OPTARG};;
s) SAMPLEID=${OPTARG};;
esac
done

source /home/${USER}/.bashrc
source activate wol

echo "PROCESSING SAMPLE: ${SAMPLEID}"
SAMPLEROOT="${ROOTDIR}$SAMPLEID"
OUTPUTROOT='/central/groups/MazmanianLab/joeB/PDMBS/classification/' 
cd $OUTPUTROOT

#_____________________________________________________________________
    
refseq_kraken_run="kraken2 --db /central/groups/MazmanianLab/joeB/Downloads/refseq_pluspf_v4/ \
--threads 2 \
--paired \
--gzip-compressed \
--classified-out "$OUTPUTROOT'RefSeqPlusPF_mapped_paired/classified_RefSeqPlusPF_'${SAMPLEID}'#.tsv'" \
--unclassified-out "$OUTPUTROOT'RefSeqPlusPF_mapped_paired/unclassified_RefSeqPlusPF_'${SAMPLEID}'#.tsv'" \
--report "$OUTPUTROOT'RefSeqPlusPF_mapped_paired/'${SAMPLEID}'ll_report_RefSeqPlusPF.tsv'" \
"${SAMPLEROOT}_1.fq.gz" \
"${SAMPLEROOT}_2.fq.gz""

echo $refseq_kraken_run
$refseq_kraken_run

#_____________________________________________________________________

uhgg_kraken_run="kraken2 --db /central/groups/MazmanianLab/joeB/Downloads/uhgg_kraken2-db/ \
--threads 2 \
--paired \
--gzip-compressed \
--classified-out "$OUTPUTROOT'UHGG_mapped_paired/classified_UHGG_'${SAMPLEID}'#.tsv'" \
--unclassified-out "$OUTPUTROOT'UHGG_mapped_paired/unclassified_UHGG_'${SAMPLEID}'#.tsv'" \
--report "$OUTPUTROOT'UHGG_mapped_paired/'${SAMPLEID}'_report_UHGG.tsv'" \
"${SAMPLEROOT}_1.fq.gz" \
"${SAMPLEROOT}_2.fq.gz""

echo $uhgg_kraken_run
$uhgg_kraken_run

#_____________________________________________________________________

WoL_kraken_run="kraken2 --db /central/groups/MazmanianLab/joeB/WebOfLife/databases/kraken2/ \
--threads 2 \
--paired \
--gzip-compressed \
--classified-out "$OUTPUTROOT'WoL_mapped_paired/classified_WoL_'${SAMPLEID}'#.tsv'" \
--unclassified-out "$OUTPUTROOT'WoL_mapped_paired/unclassified_WoL_'${SAMPLEID}'#.tsv'" \
--report "$OUTPUTROOT'WoL_mapped_paired/'${SAMPLEID}'_report_WoL.tsv'" \
"${SAMPLEROOT}_1.fq.gz" \
"${SAMPLEROOT}_2.fq.gz""

echo $WoL_kraken_run
$WoL_kraken_run