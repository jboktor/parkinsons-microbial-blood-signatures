#!/bin/bash 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=1 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=100G   # memory per CPU core
#SBATCH -J "bracken2biom"   # job name
# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=<jboktor>@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

while getopts r:s:o: option
do
case "${option}"
in
r) ROOTDIR=${OPTARG};;
s) SAMPLEID=${OPTARG};;
o) OUTPUTDIR=${OPTARG};;
esac
done

taxrank=( ["D"]="domain" ["K"]="kingdom" ["P"]="phylum" ["C"]="class" ["O"]="order" ["F"]="family" ["G"]="genus" ["S"]="species")
cd $OUTPUTROOT
kraken-biom *__species_report.tsv -o RefSeqPlusPF_bracken__species.tsv --fmt tsv