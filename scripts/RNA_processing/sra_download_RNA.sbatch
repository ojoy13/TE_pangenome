#!/bin/bash

## Example usage:
## sbatch --array 0-43 sra_download_bulk.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH --array 0-41

# Job name and output
#SBATCH -J sraDownloadBulk
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
outDir="/TE_pangenome/data/RNA/fastqs" 
sraFile="/TE_pangenome/data/RNA/sra_HGSVC3_RNA.txt" 

# Make outDir, if it doesn't exist
mkdir -p $outDir

# set query files
queries=($(cat $sraFile | grep -v "#" | awk '{print $1}'))

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Downloading sra...")

fasterq-dump -e 8 -O $outDir -v -v \
-t /Users/data/fasterq-dump-tmp \
${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")

