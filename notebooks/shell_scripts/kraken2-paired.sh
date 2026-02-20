#!/bin/bash
#______________________________________________________________________________
#                     Slurm Construction Section

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#______________________________________________________________________________

while getopts "s:r:f:" opt
do
    case "$opt" in 
        s ) SAMPLE_NAME="$OPTARG" ;;
        r ) REPORT_DIR="$OPTARG" ;;
        f ) FASTQ_DIR="$OPTARG" ;;
    esac
done

# activate enviornment
source /home/${USER}/.bashrc
source activate pdmbsR
# define enviornmental vars
threads=$SLURM_CPUS_PER_TASK 
# REPORT_DIR="/central/groups/MazmanianLab/joeB/PDMBS/workflow/WGS/results/kraken2"
# FASTQ_DIR="/central/groups/MazmanianLab/joeB/PDMBS/workflow/WGS/clean_reads"
FORWARD_RD="${FASTQ_DIR}/${SAMPLE_NAME}_R1.fastq.gz"
REVERSE_RD="${FASTQ_DIR}/${SAMPLE_NAME}_R2.fastq.gz"

date
echo "PROCESSING SAMPLE: ${SAMPLE_NAME}"
echo "FORWARD READ: ${FORWARD_RD}"
echo "REVERSE READ: ${REVERSE_RD}"
echo "Using ${threads} threads"
kraken2 --version

kraken_run() {
    echo "Reference DB: $1"
    echo "Reference DB PATH: $2"
    
    kraken2 --db $2 \
    --threads ${threads} \
    --confidence 0.3 \
    --memory-mapping \
    --paired \
    --gzip-compressed \
    --report "${REPORT_DIR}/$1_reports/${SAMPLE_NAME}_$1_report.tsv" \
    --classified-out "${REPORT_DIR}/$1_classified-reads/${SAMPLE_NAME}#_$1_classified.fastq" \
    "${FORWARD_RD}" "${REVERSE_RD}" \
    > /dev/null

    gzip "${REPORT_DIR}/$1_classified-reads/${SAMPLE_NAME}_1_$1_classified.fastq"
    gzip "${REPORT_DIR}/$1_classified-reads/${SAMPLE_NAME}_2_$1_classified.fastq"
}

# UHGG Classification
time kraken_run "UHGG" "/central/groups/MazmanianLab/joeB/Downloads/uhgg_kraken2-db/"
#"/dev/shm/uhgg_kraken2-db/"

# RefSeqPlusPF Classification
time kraken_run "RefSeqPlusPF" "/central/groups/MazmanianLab/joeB/Downloads/refseq_pluspf_v4/"
#  "/dev/shm/refseq_pluspf_v4/"
