#!/bin/bash
#______________________________________________________________________________
#                     Slurm Construction Section

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20GB
#SBATCH --time=5:00:00
#______________________________________________________________________________

while getopts "s:w:o:" opt
do
    case "$opt" in 
        s ) SAMPLE_NAME="$OPTARG" ;;
        w ) WORK_DIR="$OPTARG" ;;
        o ) OUTPUT_DIR="$OPTARG" ;;
    esac
done

# activate enviornment
source /home/${USER}/.bashrc
source activate pdmbsR
# define enviornmental vars
threads=$SLURM_CPUS_PER_TASK 
sex_fasta="/central/groups/MazmanianLab/joeB/Downloads/sex-chromosome-fragments/all_sex_sequences.fasta"
bbduk_metrics="${OUTPUT_DIR}/stats_clean_reads"
readsdir_clean="${OUTPUT_DIR}/clean_reads"
readsdir_raw="${WORK_DIR}/raw_reads"
reads_raw="${readsdir_raw}/${SAMPLE_NAME}.fastq"


# partition BAM into F/R/S fastq files
bam_to_fastq() {
    samtools bam2fq -@ ${threads} \
    "${OUTPUT_DIR}/bam_unmapped/${SAMPLE_NAME}.bam" > "${reads_raw}"
}
bam_to_fastq

# minlen=50 \ #after trimming, discard reads if this short  
# k=25 \ #kmer length  
# mink=8 \ #look for shorter kmers at read tips to this min 
# ktrim=rl \ # trim bases that problematic sequences that map to sex chromosome fragments, trim to the right
# ref="${sex_fasta}" \  # reads that align to sex-chrom. fragments
# hdist=1 \ # max hamming distance for reference kmers - num mistmatch allowed 
# qtrim=rl \ # trim R and L ends  at trimq (after kmer search) 
# trimq=20 \  # regions with average q below this will be trimmed 
# bhist="bhist.txt" \ # Base composition histogram by position. 
# qhist="qhist.txt" \ #  Quality histogram by position.
# gchist="gchist.txt" \ #  GC content
# aqhist="aqhist" \ #  Histogram of average read quality.
# lhist="lhist.txt" \ #  Read length histogram.

read_qc() {
    ~/bbmap/bbduk.sh in="${reads_raw}" \
        out="${readsdir_clean}/${SAMPLE_NAME}_R1.fastq" \
        out2="${readsdir_clean}/${SAMPLE_NAME}_R2.fastq" \
        -Xmx15g  \
        minlen=50 \
        k=25 \
        mink=8 \
        ktrim=r \
        ref="${sex_fasta}" \
        hdist=1 \
        overwrite=true \
        qtrim=rl \
        trimq=10 \
        t=${threads} \
        bhist="${bbduk_metrics}/${SAMPLE_NAME}_bhist.txt" \
        qhist="${bbduk_metrics}/${SAMPLE_NAME}_qhist.txt" \
        gchist="${bbduk_metrics}/${SAMPLE_NAME}_gchist.txt" \
        aqhist="${bbduk_metrics}/${SAMPLE_NAME}_aqhist.txt" \
        lhist="${bbduk_metrics}/${SAMPLE_NAME}_lhist.txt" \
         >& "${bbduk_metrics}/${SAMPLE_NAME}_bbduk_stdout.txt"
}
read_qc

rm ${reads_raw}
gzip "${readsdir_clean}/${SAMPLE_NAME}_R1.fastq"
gzip "${readsdir_clean}/${SAMPLE_NAME}_R2.fastq"
