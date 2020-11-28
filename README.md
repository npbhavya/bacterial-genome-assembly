# Genome assembly
work in progress 

### Description 
In this project the cultured single cell bacterial genomes were sequenced using either Illumina or Nanopore sequencing, or both. The genomes were assembled using mutliple sequencing technologies, so we applied hybrid assembly methods to assemble complete genomes. 

### Assembly 
For help with which assember to use on the dataset, this is a really useful guide https://github.com/rrwick/Trycycler/wiki/Guide-to-bacterial-genome-assembly
For the assembled genomes for this project, we had Nanopore and Illumina sequences
Assembler used - Unicycler

**Script to run assembly** 
The snakemake script runs QC and unicycler on all the genomes sequenced. Below is a description of the snakemake rules, 
- Illumina QC - Prinseq (https://github.com/Adrian-Cantu/PRINSEQ-plus-plus)
- Nanopore QC - Filtlong (https://github.com/rrwick/Filtlong)
  Note: the steps that include converting fast5 sequences to fastq was already. In the process some steps of QC was already done as well (need to make a note of the   program used here and the steps). Filtlong can be therefore redundant and can be skipped
- Assembly - Unicycler (https://github.com/rrwick/Unicycler)

### gfa graph traversal 
After assembly, the genomes can still be partially fragemented. In this case, if the number of contigs are few (worked with less than 30 for now), but the bandage plot shows that chromosomes are complete and circular, we can traverse through the gfa file to complete the genome. 

Here we applied a depth first based traversal algorithm (dfs) to walk through the graph (gfa file) that lists the contigs as nodes, and all possible connection between then as edges. The jupyter notebook contains a script that makes sure all the nodes are traveresed at least once, and the circle back to the beginning. Therefore, completing the genome. 

Further evidence from the nanopore reads is used to confirm ribosomal sequence and repeat regions.





