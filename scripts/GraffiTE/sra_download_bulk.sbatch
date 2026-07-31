#!/bin/bash

## Script for downloading HGSVC3 fastq SRA files
## Example usage:
## sbatch  sra_download_bulk_OJ.sbatch

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH --array 0-147

# Job name and output
#SBATCH -J sraDownloadBulk
#SBATCH -o /Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%A_%a.err

pwd; hostname; date

# set query files
queries=($(cat $sraFile | grep -v "#" | awk '{print $1}'))

# set variables
outDir="/TE_pangenome/GraffiTE/outputs/fastqs"
sraFile="/TE_pangenome/GraffiTE/data/sra_HGSVC3.txt"

# Make outDir, if it doesn't exist
mkdir -p $outDir

# run the fasterq-dump

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Downloading sra...")

fasterq-dump -e 8 -O $outDir -v -v \
-t /data/fasterq-tmp \
${queries[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Done!")

