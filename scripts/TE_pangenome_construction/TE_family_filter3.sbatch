#!/bin/bash

## Script for final filter keeping variants with < 3 TEs in INFO/repeat_ids
##
##
## Example usage:
## sbatch TE_family_filter3.sbatch

## General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=1:00:00
#SBATCH --mem=128GB

# Job name and output
#SBATCH -J TE_family_filter3
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
inDir="/TE_pangenome/TE_pangenome_construction/outputs/repeatmasker"
outDir="/TE_pangenome/TE_pangenome_construction/outputs/filtered_pangenome"

# Make outDir, if it doesn't exist
mkdir -p $outDir

# Filter to 3TEs per hit
echo $(date +"[%b %d %H:%M:%S] Starting...")

bcftools view $inDir/genotypes.repmasked.vcf.gz | awk '
    BEGIN {FS=OFS="\t"}
    /^#/ {print; next}  # Print all header lines and skip to next line
    {
        split($8, a, ";");
        for (i in a) {
            if (a[i] ~ /^repeat_ids=/) {
                split(a[i], b, "=");
                # Count commas in the repeat_ids value, add 1 for number of entries
                n = gsub(",", ",", b[2]);
                if (n+1 <= 3) {  # If 3 or fewer repeat IDs
                    print;
                    break;  # Break once we find and process repeat_ids
                }
            }
        }
    }
' > $outDir/pangenome3TE.vcf

echo $(date +"[%b %d %H:%M:%S] Done!")