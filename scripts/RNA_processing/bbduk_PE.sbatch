#!/bin/bash

## Script for running bbduk
##
## Example usage:
## sbatch  bbduk_PE.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --array 0-20

# Job name and output
#SBATCH -J bbduk_PE
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
inDir="/TE_pangenome/outputs/RNA/fastqs"  
outDir="/TE_pangenome/outputs/RNA/trimmed_fastqs"
adapterFile="/path/to/dir/bbmap/38.05/resources/adapters.fa"
numThreads=1

# Make outDir, if it doesn't exist
mkdir -p $outDir

# Define query files
queries=($(ls $inDir/*fastq | xargs -n 1 basename | sed 's/_R1.fastq//g' | sed 's/_R2.fastq//g' | uniq))

echo "bbduk version: "$(bbduk.sh --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Running bbduk...")

bbduk.sh -Xmx32g \
in1=$inDir/${queries[$SLURM_ARRAY_TASK_ID]}_R1.fastq \
in2=$inDir/${queries[$SLURM_ARRAY_TASK_ID]}_R2.fastq \
out1=$outDir/${queries[$SLURM_ARRAY_TASK_ID]}_R1_trimmed.fastq \
out2=$outDir/${queries[$SLURM_ARRAY_TASK_ID]}_R2_trimmed.fastq \
ref=${adapterFile} \
ktrim=r k=31 mink=11 hdist=1 tpe tbo \
qtrim=r trimq=10 \
t=${numThreads}

echo $(date +"[%b %d %H:%M:%S] Done!")
