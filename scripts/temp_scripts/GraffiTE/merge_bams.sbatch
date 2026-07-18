#!/bin/bash

## Script for running merge_bams.py to merge bams into samples
##
## Example usage:
## sbatch --array=0-20 merge_bams.sbatch

## General settings
#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --time=72:00:00
#SBATCH --mem=50GB

# Job name and output
#SBATCH -J merge_bams
#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

# Load the samtools module
module load samtools

# Define input and output directories
samples_csv="sample_ID_bam_list.csv"
bam_dir="/scratch/Users/olde5615/data/scale_up"
output_dir="/scratch/Users/olde5615/data/scale_up/bams_withID"

# Get the list of unique sample IDs
sample_ids=($(cut -d, -f1 "$samples_csv" | sort | uniq))

# Get the sample ID for this job array task
sample_id=${sample_ids[$SLURM_ARRAY_TASK_ID]}

# Run the merge for this sample
echo "Processing sample: $sample_id"
bam_files=($(grep "$sample_id" "$samples_csv" | cut -d, -f2 | sed "s|^|$bam_dir/|"))
if [ ${#bam_files[@]} -gt 1 ]; then
    merged_bam="$output_dir/$sample_id.merged.bam"
    echo "Merging ${#bam_files[@]} BAM files for $sample_id into $merged_bam..."
    samtools merge -f -o "$merged_bam" "${bam_files[@]}"
else
    echo "Only one BAM file for $sample_id. Skipping merge."
fi

echo "Done processing sample: $sample_id"