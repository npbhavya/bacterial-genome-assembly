# Genome assembly
work in progress 

### Description 
In this project the cultured single cell bacterial genomes were sequenced using either Illumina or Nanopore sequencing, or both. The genomes were assembled using mutliple sequencing technologies, so we applied hybrid assembly methods to assemble complete genomes. 

### Assembly 
For help with which assember to use on the dataset, this is a really useful guide https://github.com/rrwick/Trycycler/wiki/Guide-to-bacterial-genome-assembly

### gfa graph traversal 
After assembly, the genomes are generally still partially fragemented. In this case, if the number of contigs are few (worked with less than 30 for now), we can traverse through the gfa file to complete the genome. Here we applies a depth first based traversal algorithm (dfs) to walk through the graph (gfa file) that lists the contigs as nodes, and all possible connection between then as edges. The jupyter notebook contains a script that makes sure all the nodes are traveresed at least once, and the circle back to the beginning. Therefore, completing the genome. 

Further evidence from the nanopore reads is used to confirm ribosomal sequence and repeat regions. 





