#!/bin/bash

#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --time=72:00:00
#SBATCH --mem=16GB
#SBATCH -J fc_perSample
#SBATCH --array=1-21%6

#SBATCH -o /scratch/Users/%u/data/slurmOut/slurm-%A_%a.out
#SBATCH -e /scratch/Users/%u/data/slurmErr/slurm-%A_%a.err

pwd; hostname; date

module load subread

inDir=/scratch/Users/olde5615/data/graph21_RNA
outDir=${inDir}/featureCounts_transcripts_gene_id/perSample
gtf=/scratch/Users/olde5615/data/gencode.v38.ComprehensiveAnnotation.Main.gtf
feature=exon
strandOption=2
bamlist=${inDir}/featureCounts_transcripts/bamlist.txt


mkdir -p "$outDir"


bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$bamlist")
base=$(basename "$bam" .bam)

echo "$(date +"[%b %d %H:%M:%S]") Starting"


echo "Running featureCounts on $bam"
featureCounts --verbose -T "$SLURM_CPUS_PER_TASK" -O -p -s "$strandOption" \
  -t "$feature" -g gene_id \
  -a "$gtf" \
  -o "$outDir/${base}_${strandOption}.featureCounts.txt" \
  "$bam"
 
 echo "$(date +"[%b %d %H:%M:%S]") Done"
