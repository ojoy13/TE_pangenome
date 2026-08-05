## Background:
Transposon Element (TE) focused pangenome pipeline for polymorphic transposable elements in diverse populations.
Specifically: Alu, LINE1, SVA, & HERV
Matched short-read RNA is then analyzed using featureCounts and DESeq2. The alignment is TE-pangenome variant-informed and looks for nearby genes that may also be impacted by gene dysregulation.

Look out for publication through microPublication!
(Coming soon)

### Abstract:
Using a pangenome allows us to identify uncommon transposable elements (TEs) that are difficult to find with current alignment methods. We know that TEs are one source of interindividual variation, but we are still unsure of the impact of TE variation on phenotypes. Here we use a TE-focused pangenome approach to find TEs associated with gene dysregulation in 21 diverse genomes. We found that TE polymorphisms can sometimes overlap with introns and enhancers, providing more insight into why they lead to differential phenotype outcomes. 

##  Programs used:
#### TE pangenome construction:
- fasterqdump v3.0.0 (https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
- bbmap (bbDuk) v38.05
- minimap2 v2.17-r941 (https://github.com/lh3/minimap2)
- NanoStat v1.6.0   (https://pypi.org/project/NanoStat/)
    - python v3.10.8
- NanoPlot v1.39.0 (optional) (https://github.com/wdecoster/NanoPlot)
    - python v3.10.8
- samtools v1.16.1 (https://github.com/samtools/samtools/releases/)
    - htslib v1.16
- bcftools v1.8 (https://github.com/samtools/bcftools/releases/)
- sniffles v2.4 (https://github.com/fritzsedlazeck/Sniffles)
    - python v3.10.8 (https://www.python.org/downloads/release/python-3108/)
- RepeatMasker v4.1.0 (https://www.repeatmasker.org/RepeatMasker/)
- R v4.4 https://cran.r-project.org/bin/windows/base/old/4.4.0/
- bcftools v1.8 https://github.com/samtools/bcftools/releases/tag/1.8
- bedtools v2.28.0  https://github.com/arq5x/bedtools2/releases/tag/v2.28.0
- samtools v1.16    https://github.com/samtools/samtools/releases/tag/1.16

### RNA alignment:
- fasterqdump v3.0.0 (https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
- bbmap (bbDuk) v38.05
- NanoStat v1.6.0   (https://pypi.org/project/NanoStat/)
    - python v3.10.8 (https://www.python.org/downloads/release/python-3108/)
- NanoPlot v1.39.0 (optional) (https://github.com/wdecoster/NanoPlot)
    - python v3.10.8 (https://www.python.org/downloads/release/python-3108/)
- minimap2 v2.17-r941 (https://github.com/lh3/minimap2)

### FeatureCounts:
- Hisat2 v2.1.0
- bamCoverage v3.0.1
- featureCounts v1.6.2

### DESeq2:
- R 4.4 https://cran.r-project.org/bin/windows/base/old/4.4.0/
- BiocManager v3.19
    - DESeq2 v1.44.0

## Workflow
See /TE_pangenome/workflow/ for full diagrams of each step!
