#!/bin/bash

## Script for identifying structural variants by POPULATION call
## Note: Sniffles2 requires python 3
## Example usage:
## sbatch sniffles_popn.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=1:00:00
#SBATCH --mem=256GB

# Job name and output
#SBATCH -J sniffles_popn_call
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
inDir="/TE_pangenome/TE_pangenome_construction/outputs/sniffles" 
outDir="/TE_pangenome/TE_pangenome_construction/outputs/sniffles_popn"
ref="/path/to/dir/hg38.main" 

# Make outDir, if it doesn't exist
mkdir -p $outDir

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Concatenate snf files to tsv...")

ls $inDir/*.snf > $inDir/snfs.tsv

echo $(date +"[%b %d %H:%M:%S] Starting sniffles for population SVs...")

sniffles --minsvlen 100  --threads 16 --reference $ref --input $inDir/snfs.tsv --vcf $outDir/genotypes_unfiltered.vcf

echo $(date +"[%b %d %H:%M:%S] Starting bcftools filtering...")
  
bcftools filter --threads 16 -i 'INFO/SVTYPE == "INS" | INFO/SVTYPE == "DEL"' $outDir/genotypes_unfiltered.vcf | awk '$5 != "<INS>" && $5 != "<DEL>"' > $outDir/sniffles2_variants.vcf

echo $(date +"[%b %d %H:%M:%S] Done!")
