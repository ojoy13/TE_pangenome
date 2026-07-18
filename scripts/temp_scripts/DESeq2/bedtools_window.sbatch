#!/bin/bash

## Script looking at neighboring genes to TE Polymorphisms
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
#SBATCH -J bed_win
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%j.err

module load bedtools

echo $(date +"[%b %d %H:%M:%S]") "Done!"

bedtools window -w 10000 -a /scratch/Users/olde5615/data/pangenome21_phased/21_graphs_31MAY26/phased_vcf/pangenome_21_phased_nosymbolic.vcf.gz -b gtf_transcripts_bed_output.bed > gtf_transcripts_bedwin_10kb.bed

echo $(date +"[%b %d %H:%M:%S]") "Done!"
