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
- minimap2 v2.17-r941 (https://github.com/lh3/minimap2)
- samtools v1.16.1 (https://github.com/samtools/samtools/releases/)
    - htslib v1.16
- sniffles v2.4 (https://github.com/fritzsedlazeck/Sniffles)
    - python v3.10.8 (https://www.python.org/downloads/release/python-3108/)
- NanoPlot v1.39.0 (optional) (https://github.com/wdecoster/NanoPlot)
    - python v3.10.8
- RepeatMasker v4.1.0 (https://www.repeatmasker.org/RepeatMasker/)
- bcftools v1.8 (https://github.com/samtools/bcftools/releases/)

### RNA alignment:
-
-
### FeatureCounts:
-
-
-
### DESeq2:
-
-
-
## Workflow
See /TE_pangenome/workflow/ for full diagrams of each step!
