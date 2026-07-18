#!/bin/bash

## Script for identifying structural variants by POPULATION call
##
## Example usage:
## inDir=bams_withID outDir=. sbatch sniffles_popn.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=1:00:00
#SBATCH --mem=256GB

# Job name and output
#SBATCH -J sniffles_popn_call
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%j.err

#activate virtual environment that contains sniffles
source ~/jhub_venv_3.10.8/bin/activate

#modules
module load bcftools
module load samtools

# run the thing
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Concatenate snf files to tsv...")

echo $(date +"[%b %d %H:%M:%S] Starting sniffles for population SVs...")

sniffles --minsvlen 100  --threads 16 --reference /Shares/CL_Shared/db/genomes/hg38/fa/hg38.main.fa --input $inDir/snfs.tsv --vcf $outDir/genotypes_unfiltered.vcf

echo $(date +"[%b %d %H:%M:%S] Starting bcftools filtering...")
  
bcftools filter --threads 16 -i 'INFO/SVTYPE == "INS" | INFO/SVTYPE == "DEL"' $inDir/genotypes_unfiltered.vcf | awk '$5 != "<INS>" && $5 != "<DEL>"' > $outDir/sniffles2_variants.vcf

echo $(date +"[%b %d %H:%M:%S] Done!")
