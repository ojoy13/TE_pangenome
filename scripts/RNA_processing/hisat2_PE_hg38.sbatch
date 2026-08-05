#!/bin/bash

## Script for running hisat2
##
## Example usage:
## sbatch  hisat2_PE.sbatch

# important note: 
# for single-ended data, rnaStrandness is R (for reverse) or F (for forward)
# for paired-ended data, rnaStrandness is RF (for reverse) or FR (for forward)

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=24:00:00
#SBATCH --mem=500GB
#SBATCH --array 0-20

# Job name and output
#SBATCH -J hisat2_PE
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
inDir="/TE_pangenome/outputs/RNA/bams_merged" 
outDir="/TE_pangenome/outputs/RNA/hisat2" 
hisatIdxDir="/path/to/dir/hisat2" 
hisatIdx="/path/to/dir/hg38.main" 
rnaStrandness=RF
numThreads=16

# Make outDir, if it doesn't exist
mkdir -p $outDir

# Set query files
# note: remove R1/R2 to generate a unique identifier for each pair of files
queries=($(ls $inDir/*_trimmed.fastq | xargs -n 1 basename | sed 's/_R1_trimmed.fastq//g' | sed 's/_R2_trimmed.fastq//g' | uniq))

echo "Starting HISAT2 alignment..."
echo "HISAT2 version: "$(hisat2 --version)
echo "samtools version: "$(samtools --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}

hisat2 --no-softclip \
--rna-strandness ${rnaStrandness} \
-p ${numThreads} \
-x ${hisatIdxDir}/${hisatIdx} \
-1 $inDir/${queries[$SLURM_ARRAY_TASK_ID]}_R1_trimmed.fastq \
-2 $inDir/${queries[$SLURM_ARRAY_TASK_ID]}_R2_trimmed.fastq \
| samtools view -q 10 -Sb - \
| samtools sort -@ ${numThreads} - -o $outDir/${queries[$SLURM_ARRAY_TASK_ID]}.bam

# Create a bai index file for upload to UCSC browser
samtools index $outDir/${queries[$SLURM_ARRAY_TASK_ID]}.bam

echo $(date +"[%b %d %H:%M:%S] Done!")
