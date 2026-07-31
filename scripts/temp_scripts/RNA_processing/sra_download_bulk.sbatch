#!/bin/bash

## Example usage:
## sbatch --array 0-43 sra_download_bulk.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J sraDownloadBulk
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
outDir=. 
sraFile=sra.txt 

# set query files
queries=($(cat $sraFile | grep -v "#" | awk '{print $1}'))

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Downloading sra...")

# don't module load sra
# since the newest version (sra3.0.0) is in my ~/.bash_profile PATH

fasterq-dump -e 8 -O $outDir -v -v \
-t /scratch/Users/ativ2716/data/fasterq-dump-temp \
${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")

