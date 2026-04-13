#!/bin/bash
#______________________________________________________________________________
#                     Slurm Construction Section

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=30:00
#______________________________________________________________________________

while getopts "s:i:o:" opt
do
    case "$opt" in 
        s ) SAMPLE_NAME="$OPTARG" ;;
        i ) INPUT_DIR="$OPTARG" ;;
        o ) OUTPUT_DIR="$OPTARG" ;;
    esac
done

# activate enviornment
source /home/${USER}/.bashrc
source activate pdmbsR
# define enviornmental vars
# INPUT_DIR="/resnick/scratch/jbok/PDMBS/WGS/BAM-unmapped-raw"
# OUTPUT_DIR="/resnick/groups/MazmanianLab/jboktor/PDMBS/workflow/WGS"
bam_filtered="${OUTPUT_DIR}/bam_unmapped/${SAMPLE_NAME}.bam"
flagstat_dir="${OUTPUT_DIR}/flagstats/03_BAM-unmapped-3328-filtered"
sequence_meta_dir="${OUTPUT_DIR}/sequence_metadata"
threads=$SLURM_CPUS_PER_TASK 

# filter PCR duplicates (1024); secondary (256) and supplementary aligned reads (2048)
samtools_filter() {
    samtools view -@ "${threads}" -h -b -F 3328 \
    "${INPUT_DIR}/${SAMPLE_NAME}_unmapped.bam" > "${bam_filtered}"
}
samtools_filter

# Collect flagstat metrics after filtering unwanted reads
flagstat_cmd() {
    samtools flagstat -@ "${threads}" -O tsv \
    "${bam_filtered}" \
    > "${flagstat_dir}/${SAMPLE_NAME}_unmapped-F3328_flagstat.tsv"
}
flagstat_cmd

# Collect sample metadata
sequence_meta() {
    samtools view -@ "${threads}" -H "${bam_filtered}" | \
    awk '$1=="@RG" {print $0}' \
    > "${sequence_meta_dir}/${SAMPLE_NAME}_sequencing_params.tsv"
}
sequence_meta
