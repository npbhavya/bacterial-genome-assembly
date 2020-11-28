# Genome assembly
work in progress 

### Description 
In this project the cultured single cell bacterial genomes were sequenced using either Illumina or Nanopore sequencing, or both. The genomes were assembled using mutliple sequencing technologies, so we applied hybrid assembly methods to assemble complete genomes. 

### Dependecies 
  - miniconda https://docs.conda.io/en/latest/miniconda.html
  - Jupyter notebook  https://jupyter.org/
  - Bandage plots https://rrwick.github.io/Bandage/   
  - within conda env 
    - unicycler  https://anaconda.org/bioconda/unicycler
    - prinseq++  https://anaconda.org/bioconda/prinseq-plus-plus
    - filtlong   https://anaconda.org/bioconda/filtlong
    - quast      https://anaconda.org/bioconda/quast
    - snakemake  https://anaconda.org/bioconda/snakemake
   
### Input files 
The input files should be placed in two directories 
  - ED-Illumina-MIGS - place the Illumina forward and reverse reads in this directory
  - ED-ONT - place all the Nanopore reads in fastq format here
  
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

Run the script by 

    `souce activate <conda env>`
  
    `snakemake -s snakefile-ED`
  
  To run a specific rule within snakefile
  
    `snakemake -s snakefile-ED <rulename>`

The above steps were run on a cluster. After assembly, the assembly.gfa plots were visualized on Banadge plots to confirm if the assembly formed complete, circular genomes. 

Assembly statictics can be run uaing quast,

    `quast.py assembly.fasta -o <assembly-quast-output>` 

### gfa graph traversal 
After assembly, the genomes can still be partially fragemented. In this case, if the number of contigs are few (worked with less than 30 for now), but the bandage plot shows that chromosomes are complete and circular, we can traverse through the gfa file to complete the genome. 

Here we applied a depth first based traversal algorithm (dfs) to walk through the graph (gfa file) that lists the contigs as nodes, and all possible connection between then as edges. The jupyter notebook contains a script that makes sure all the nodes are traveresed at least once, and the circle back to the beginning. Therefore, completing the genome. 

**Run gfa traversal**

    `jupyter notebook`
  
Change the name of the chromosome running through the python script and the output name. 
Run the script for only one chromosome within a genome at a time. If the genome includes multiple chromosomes they need to be separated out different fasta files and each run through the script separately. 
The final circular chromosomes need to be concatenated to one file- FINAL assembly

### Updates to add to the scripts
- hardcoded directory names in the snakemake file
- add quast step to the snakemake file, and the adding all the different genome assemblies quast report to one output
- add coverage calculation from illumina and nanopore assemblies




