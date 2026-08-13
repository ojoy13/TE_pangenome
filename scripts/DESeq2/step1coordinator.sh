#!/bin/sh
Rscript /scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts/deseq_init_master_sniffles.R

CSV_FILE="/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts/mini_candidate_sniffles35.csv"

# Skip the header and loop through each line of the CSV file
tail -n +2 "$CSV_FILE" | while IFS=',' read -r index ensembl chr symbol
do
    # Remove any leading or trailing spaces (if any)
    symbol=$(echo "$symbol" | xargs)

    # Export SYMBOL as a variable to pass to the sbatch script
    export a_sniffle="$symbol"

    # Submit the job via sbatch, passing the variable agenename
    sbatch --export=a_sniffle "/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts/step1_gene_Submit_Rscript.sbatch"

    # Optional: Print the symbol being passed
    echo "Submitting job for a_sniffle: $symbol"
#break #once the whole thing works for the first gene--- remove this break. It only lets the script loop once until then. 
done

# Next merge jobs

