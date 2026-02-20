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

while getopts "s:r:f:c:" opt
do
    case "$opt" in 
        s ) SAMPLE_NAME="$OPTARG" ;;
        r ) REPORT_DIR="$OPTARG" ;;
        f ) FASTQ_DIR="$OPTARG" ;;
        c ) CLASS_FQ_DIR="$OPTARG" ;;
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
    
    krakenuniq --db "/central/groups/MazmanianLab/joeB/Downloads/RefDBs/KrakenUniq/MicrobialDB" \
    --threads ${threads} \
    --preload \
    --paired \
    --output "${REPORT_DIR}/${SAMPLE_NAME}_KrakenUniq_report.tsv" \
    --only-classified-output \
    --classified-out "${CLASS_FQ_DIR}/${SAMPLE_NAME}_classified.fastq" \
    "${FORWARD_RD}" "${REVERSE_RD}" \
    > /dev/null

    gzip "${CLASS_FQ_DIR}/${SAMPLE_NAME}_classified.fastq"
}

# KrakenUniq Classification
time kraken_run
