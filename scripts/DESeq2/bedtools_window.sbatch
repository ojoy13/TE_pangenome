#!/bin/bash

## Script looking at neighboring genes to TE Polymorphisms
## gtf_geneID_bed_output.bed created by getting enhancers from gencode.v38
##
## Example usage:
## bedtools_window.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=8:00:00
#SBATCH --mem=8GB

# Job name and output
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
inDir="/TE_pangenome/outputs/TE_pangenome_construction/filtered_pangenome"
outDir="/TE_pangenome/outputs/DESEq2"

echo $(date +"[%b %d %H:%M:%S]") "Starting..."

# Run bedtools window
bedtools window -w 10000 -a $inDir/TE_pangenome_21samples.vcf.gz -b gtf_geneID_bed_output.bed > $outDir/gtf_genes_bedwin_10kb.bed

echo $(date +"[%b %d %H:%M:%S]") "Done!"