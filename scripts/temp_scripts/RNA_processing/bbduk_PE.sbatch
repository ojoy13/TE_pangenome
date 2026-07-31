#!/bin/bash

## Script for running bbduk
##
## Example usage:
## inDir=fastq outDir=fastq sbatch --array 0-21 bbduk_PE.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --mem=128G

# Job name and output
#SBATCH -J bbduk_PE
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

# Set constant variables
numThreads=1

# Load module
module load bbmap

# Define query files
# Removes R1/R2 to generate a unique identifier for each pair of files
queries=($(ls $inDir/*fastq | xargs -n 1 basename | sed 's/_R1.fastq//g' | sed 's/_R2.fastq//g' | uniq))

# define key variables
adapterFile=/opt/bbmap/38.05/resources/adapters.fa

# Run bbduk
pwd; hostname; date

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
