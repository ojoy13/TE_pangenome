import os
import csv
import subprocess
from collections import defaultdict

# Path to the sample_ID_bam_list.csv file: sample_ID, sorted.bam file
samples_csv = "/TE_pangenome/data/TE_pangenome_construction/sample_ID_bam_list.csv"

# Directory containing sorted BAM files
bam_dir = "/TE_pangenome/outputs/TE_pangenome_construction/bams_bySRA/"  

# Output directory for merged BAM files
output_dir = "/TE_pangenome/outputs/TE_pangenome_construction/bams_merged/"
os.makedirs(output_dir, exist_ok=True)

# Dictionary to group BAM files by Sample_ID
sample_to_bams = defaultdict(list)

# Read the samples.csv file
with open(samples_csv, "r") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        sample_id = row["Sample_ID"]
        bam_file = row["BAM_Files"]
        sample_to_bams[sample_id].append(os.path.join(bam_dir, bam_file))

# Merge BAM files for each Sample_ID
for sample_id, bam_files in sample_to_bams.items():
    if len(bam_files) > 1:  
        merged_bam = os.path.join(output_dir, f"{sample_id}.merged.bam")
        print(f"Merging {len(bam_files)} BAM files for {sample_id} into {merged_bam}...")
        subprocess.run(["samtools", "merge", "-o", merged_bam] + bam_files, check=True)
    else:
        print(f"Only one BAM file for {sample_id}. Skipping merge.")



