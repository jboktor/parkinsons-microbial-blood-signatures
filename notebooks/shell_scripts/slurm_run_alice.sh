#!/bin/bash
#SBATCH --job-name=run_ALICE
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=600GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=jboktor@caltech.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=/central/groups/MazmanianLab/joeB/PDMBS/parkinsons-microbial-blood-signatures/.cluster_runs/ALICE_tcrb_%j.log

source /home/${USER}/.bashrc
source activate pdmbsR

# Set paths
SCRIPT="/central/groups/MazmanianLab/joeB/software/ALICE/execute_PD_ALICE.R"

cd /central/groups/MazmanianLab/joeB/software/ALICE/

echo "Starting R script..."
Rscript "$SCRIPT"
status=$?

if [[ $status -eq 0 ]]; then
  echo "✅ Script completed successfully at $(date)."
else
  echo "❌ Script failed with status $status at $(date)."
fi
