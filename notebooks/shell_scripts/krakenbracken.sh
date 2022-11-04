#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=3:00:00   # walltime
#SBATCH --ntasks=2 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=50G   # memory per CPU core
#SBATCH -J "Krak2"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=<jboktor>@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output=/central/scratch/jbok/slurmdump/%j.out


while getopts r:s: option
do
case "${option}"
in
r) ROOTDIR=${OPTARG};;
s) SAMPLEID=${OPTARG};;
esac
done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

source /home/${USER}/.bashrc
source activate wol

echo "PROCESSING SAMPLE: "$SAMPLEID
SAMPLEROOT="${ROOTDIR}$SAMPLEID"
OUTPUTROOT='/central/groups/MazmanianLab/joeB/PDBM/classification/' 
declare -A taxrank
taxrank=( ["D"]="domain" ["K"]="kingdom" ["P"]="phylum" ["C"]="class" ["O"]="order" ["F"]="family" ["G"]="genus" ["S"]="species")
cd $OUTPUTROOT

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if [ ! -f $OUTPUTROOT'RefSeqPlusPF_mapped/'${SAMPLEID}'__report_RefSeqPlusPF.tsv' ]
then
    echo ${SAMPLEID}"__report_RefSeqPlusPF.tsv does not exist .... running script"
    
refseq_kraken_run="kraken2 --db /central/groups/MazmanianLab/joeB/Downloads/refseq_pluspf_v4/ \
--threads 2 \
--gzip-compressed \
--classified-out "$OUTPUTROOT'RefSeqPlusPF_mapped/'${SAMPLEID}'__classified_RefSeqPlusPF.tsv'" \
--unclassified-out "$OUTPUTROOT'RefSeqPlusPF_mapped/'${SAMPLEID}'__unclassified_RefSeqPlusPF.tsv'" \
--report "$OUTPUTROOT'RefSeqPlusPF_mapped/'${SAMPLEID}'__report_RefSeqPlusPF.tsv'" \
${SAMPLEROOT}".fq.gz""
echo $refseq_kraken_run
$refseq_kraken_run

# mkdir -p RefSeqPlusPF_mapped/bracken
# for rank in "${!taxrank[@]}"
# do 
# echo "$rank - ${taxrank[$rank]}"
# bracken -d /central/groups/MazmanianLab/joeB/Downloads/refseq_pluspf_v4/ \
# -i $OUTPUTROOT'RefSeqPlusPF_mapped/'${SAMPLEID}'__report_RefSeqPlusPF.tsv' \
# -o $OUTPUTROOT'RefSeqPlusPF_mapped/bracken/'${SAMPLEID}'__'${taxrank[$rank]}'.tsv' \
# -w $OUTPUTROOT'RefSeqPlusPF_mapped/bracken/'${SAMPLEID}'__'${taxrank[$rank]}'_report.tsv' \
# -l $rank -t 0
# done


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

uhgg_kraken_run="kraken2 --db /central/groups/MazmanianLab/joeB/Downloads/uhgg_kraken2-db/ \
--threads 2 \
--gzip-compressed \
--classified-out "$OUTPUTROOT'UHGG_mapped/'${SAMPLEID}'__classified_UHGG.tsv'" \
--unclassified-out "$OUTPUTROOT'UHGG_mapped/'${SAMPLEID}'__unclassified_UHGG.tsv'" \
--report "$OUTPUTROOT'UHGG_mapped/'${SAMPLEID}'__report_UHGG.tsv'" \
${SAMPLEROOT}".fq.gz""
echo $uhgg_kraken_run
$uhgg_kraken_run

# mkdir -p UHGG_mapped/bracken
# for rank in "${!taxrank[@]}"
# do 
# echo "$rank - ${taxrank[$rank]}"
# bracken -d /central/groups/MazmanianLab/joeB/Downloads/uhgg_kraken2-db/ \
# -i $OUTPUTROOT'UHGG_mapped/'${SAMPLEID}'__report_UHGG.tsv' \
# -o $OUTPUTROOT'UHGG_mapped/bracken/'${SAMPLEID}'__'${taxrank[$rank]}'.tsv' \
# -w $OUTPUTROOT'UHGG_mapped/bracken/'${SAMPLEID}'__'${taxrank[$rank]}'_report.tsv' \
# -l $rank -t 0
# done


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

WoL_kraken_run="kraken2 --db /central/groups/MazmanianLab/joeB/WebOfLife/databases/kraken2/ \
--threads 2 \
--gzip-compressed \
--classified-out "$OUTPUTROOT'WoL_mapped/'${SAMPLEID}'__classified_WoL.tsv'" \
--unclassified-out "$OUTPUTROOT'WoL_mapped/'${SAMPLEID}'__unclassified_WoL.tsv'" \
--report "$OUTPUTROOT'WoL_mapped/'${SAMPLEID}'__report_WoL.tsv'" \
${SAMPLEROOT}".fq.gz""
echo $WoL_kraken_run
$WoL_kraken_run

# mkdir -p WoL_mapped/bracken
# for rank in "${!taxrank[@]}"
# do 
# echo "$rank - ${taxrank[$rank]}"
# bracken -d /central/groups/MazmanianLab/joeB/WebOfLife/databases/bracken/ \
# -i $OUTPUTROOT'WoL_mapped/'${SAMPLEID}'__report_WoL.tsv' \
# -o $OUTPUTROOT'WoL_mapped/bracken/'${SAMPLEID}'__'${taxrank[$rank]}'.tsv' \
# -w $OUTPUTROOT'WoL_mapped/bracken/'${SAMPLEID}'__'${taxrank[$rank]}'_report.tsv' \
# -l $rank -t 0
# done


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
else
    echo "Found: "$OUTPUTROOT'RefSeqPlusPF_mapped/'${SAMPLEID}'__report_RefSeqPlusPF.tsv'
fi




