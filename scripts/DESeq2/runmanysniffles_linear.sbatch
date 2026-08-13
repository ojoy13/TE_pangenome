#!/bin/bash
## Requires: Rv4.4.0

#SBATCH --job-name=snif_transcript
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=01:00:00
#SBATCH --partition=short
#SBATCH --mem=3GB

#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err
#SBATCH --array=1-3600%32   # <-- TEMP: will override on command line

pwd; hostname; date

rscript_path=/TE_pangenome/scripts/DESEq2/deseq_1sniffle_forlooping.R
csv_file=/TE_pangenome/outputs/DESEq2/just_candidate_sniffles.csv

printf "\nR code being run: %s" "${rscript_path}"
printf "\nRun on: %s" "$(hostname)"
printf "\nRun from: %s" "$(pwd)"
printf "\nScript: %s\n" "$0"
date
printf "\nYou have requested %s core(s).\n" "$SLURM_CPUS_ON_NODE"
printf "\nSLURM_ARRAY_TASK_ID = %s\n" "$SLURM_ARRAY_TASK_ID"

# Map array index to the correct CSV line (skip header)
line_num=$((SLURM_ARRAY_TASK_ID + 1))

# Extract that line from the CSV
line=$(sed -n "${line_num}p" "$csv_file")

# If the file is shorter than the array range, line will be empty
if [[ -z "$line" ]]; then
    printf "\nNo data line for index %s (line %s in %s)\n" \
           "$SLURM_ARRAY_TASK_ID" "$line_num" "$csv_file"
    exit 0
fi

# If there are multiple columns, split them; here we take the first as a_sniffle
IFS=, read -r a_sniffle _ <<< "$line"

printf "\nRunning sniffle: %s\n" "$a_sniffle"

Rscript "${rscript_path}" $a_sniffle

echo "Job finished at $(date +"%T %a %d %b %Y")"
