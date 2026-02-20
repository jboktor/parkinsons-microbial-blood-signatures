#!/bin/bash

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                     Slurm Construction Section

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL

while getopts "i:o:" opt
do
    case "$opt" in
        i ) URL="$OPTARG" ;;
        o ) OUTPUTDIR="$OPTARG" ;;
    esac
done

# activate enviornment
source /home/${USER}/.bashrc
source activate pdmbsR
echo "$URL"
echo "$OUTPUTDIR"
mkdir -p  $OUTPUTDIR

# command to download file from google cloud
download_cmd="gsutil cp "${URL}" "${OUTPUTDIR}""
echo $download_cmd
$download_cmd

