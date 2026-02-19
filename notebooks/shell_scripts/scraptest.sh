#!/bin/bash

while getopts "i:o:" opt
do
    case "$opt" in
        i ) URL="$OPTARG" ;;
        o ) OUTPUTDIR="$OPTARG" ;;
    esac
done

source /home/${USER}/.bashrc
source activate pdmbsR

# Begin script in case all parameters are correct
echo "$URL"
echo "$OUTPUTDIR"

# command to download file from google cloud
download_cmd="gsutil cp "${URL}" "${OUTPUTDIR}""
echo $download_cmd
gsutil cp "$URL" "$OUTPUTDIR"
