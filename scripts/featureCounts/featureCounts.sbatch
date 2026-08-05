#!/bin/bash

#SBATCH -p long
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --time=72:00:00
#SBATCH --mem=16GB
#SBATCH -J fc_perSample
#SBATCH --array=1-21%6

#SBATCH -J featureCounts
#SBATCH -o /Users/%u/data/slurmOut/slurm-%j.out
#SBATCH -e /Users/%u/data/slurmErr/slurm-%j.err

pwd; hostname; date

# Set variables
inDir="/TE_pangenome/outputs/RNA/bigWigs"
outDir="/TE_pangenome/outputs/featureCounts/perSample"
gtf="/path/to/dir/gencode.v38.ComprehensiveAnnotation.Main.gtf"
feature=exon
strandOption=2
bamlist="/TE_pangenome/data/featureCounts/bamlist_sorted.txt"


mkdir -p "$outDir"


bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$bamlist")
base=$(basename "$bam" .bam)

echo "$(date +"[%b %d %H:%M:%S]") Starting"


echo "Running featureCounts on $bam"
featureCounts --verbose -T "$SLURM_CPUS_PER_TASK" -O -p -s "$strandOption" \
  -t "$feature" -g transcript_id \
  -a "$gtf" \
  -o "$outDir/${base}_${strandOption}.featureCounts.txt" \
  "$bam"
 
 echo "$(date +"[%b %d %H:%M:%S]") Done"
