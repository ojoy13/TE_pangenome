#!/bin/bash

## Example usage:
## outDir=. sraFile=sra_DNA.txt sbatch --array 0-4 sra_download.sbatch

# SBATCH settings (unchanged)
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=72:00:00
#SBATCH --mem=5GB
#SBATCH -J sraDownloadBulk
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

# Load modules (if needed)
#module load sratoolkit/3.0.0  # or your version

# Set query files and dbGaP key
queries=($(cat $sraFile | grep -v "#" | awk '{print $1}'))
dbgap_key=/Users/olde5615/dbGaP_keys/prj_31282_D38373.ngc

# Debug info
pwd; hostname; date
echo "Processing file: ${queries[$SLURM_ARRAY_TASK_ID]}"
echo "Using dbGaP key: $dbgap_key"


# Download SRA
echo "$(date +"[%b %d %H:%M:%S] Converting to FASTQ...")"
fasterq-dump -e 8 --ngc "$dbgap_key" -O "$outDir" -t /scratch/Users/olde5615/data/fasterq-dump-temp \
  --split-files "${queries[$SLURM_ARRAY_TASK_ID]}"


echo "$(date +"[%b %d %H:%M:%S] Done!")
