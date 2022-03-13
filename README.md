[![DOI](https://zenodo.org/badge/291788219.svg)](https://zenodo.org/badge/latestdoi/291788219)

# Genome assembly

### Description 
In this project the cultured single cell bacterial genomes were sequenced on both Illumina and Nanopore sequencing. This repository holds the scripts and commands used to assemble these bacterial genomes to complete circualr chromosomes. The scripts are saved as job scripts that can be submitted to SLURM job scheduler on high performace clusters and as snakemake workflow. 

### Installation

#### Dependecies 
  - miniconda https://docs.conda.io/en/latest/miniconda.html
  - Jupyter notebook  https://jupyter.org/
  - Bandage plots https://rrwick.github.io/Bandage/   
  - within conda env 
  `conda install -c bioconda prinseq` \
  `conda install -c bioconda unicycler` \
  `conda install -c bioconda flye` \
  `conda install -c bioconda filtlong` \
  `conda install -c conda-forge mamba` \
  `conda install -c bioconda snakemake` \
  `conda install -c conda-forge raven` \
  `conda install -c bioconda miniasm` \
  `pip install quast`
   
### Input files 
The input files should be placed in two directories 
  - Illumina - place the Illumina forward and reverse reads in this directory
  - ONT - place all the Nanopore reads in fastq format here

### Genome assembly workflow
All the commands are not yet automated. For now the breakdown of the steps is shown below 
- **Guppy basecaller** to convert fast5 sequences to fastq sequences - saved as a script "guppy_basecaller.sh"
- **snakemake workflow** that runs QC and assembly using unicycler - saved to "snakemake-ED", and "QC-And-Assembly.sh" script
- **Visualize the graph** file generated from unicycler on Bandage application
- Try **other assemblers** in case the genomes are fragmented - script saved to "snakefile-other-assemblers", and "other-assembler.sh" script
- **Graph traversal script** to build a complete geneome - saved as jupyter notebook "graph_traversal.ipynb"
- **Polishing the assembly** polish the complete genomes with medaka and Pilon - scripts saved to "medaka.sh", and "pilon.sh"
- **Assembly quality assessment** - run quast and ideel - saved to script "ideel.sh"

The job scripts have to be filled out with the filepaths and any changes to the folder names as neceesary 

### Assembly 
For help with which assembler to use on the dataset, this is a really useful guide https://github.com/rrwick/Trycycler/wiki/Guide-to-bacterial-genome-assembly
For the assembled genomes for this project, we had both nanopore and illumina sequences with relatively high coverage of the genomes
Primary assembler used - Unicycler

**Script to run assembly** 
The snakemake script runs QC first 
- Illumina QC - Prinseq (https://github.com/Adrian-Cantu/PRINSEQ-plus-plus)
- Nanopore QC - Filtlong (https://github.com/rrwick/Filtlong)
- Assembly - Unicycler (https://github.com/rrwick/Unicycler)

Run the script by 

    souce activate genome assembly
    snakemake -s snakefile-ED
  
  To run a specific rule within snakefile
  
    snakemake -s snakefile-ED <rulename>

Assembly statictics can be run using quast,

    quast.py assembly.fasta -o <assembly-quast-output> 

After assembly, the assembly.gfa plots were visualized on Banadge plots to confirm if the assembly formed complete, circular genomes. Below is an example of a good assembly, where the graph shows one circular contig was assembled.  
![image](https://user-images.githubusercontent.com/8825721/158039332-fddca173-8fa5-427f-97bc-b851ba8942cb.png)


### gfa graph traversal 
In most cases, the genomes are still partially fragmented. In cases where the fragmets looks simple to resolve, a python script was written to travese the path and build a complete genome assembly. The script uses a depth first based traversal algorithm (dfs) to walk through the graph (gfa file) that lists the contigs as nodes, and all possible connection between then as edges. The jupyter notebook contains a script that makes sure all the nodes are traveresed at least once, and the circle back to the beginning.

<img width="313" alt="image" src="https://user-images.githubusercontent.com/8825721/158039558-e0116547-18c4-4448-8a66-49c4558d46e3.png">

**Run gfa traversal script on notebook**

    jupyter notebook
  
Change the name of the chromosome running through the python script and the output name. 
Run the script for only one chromosome within a genome at a time. If the genome includes multiple chromosomes they need to be separated out different fasta files and each run through the script separately. 
The final circular chromosomes need to be concatenated to one file- FINAL assembly

### Polishing the assembly 
The final assembled genome is polished with evidence from both the nanopore and illumina reads. 
**Polishing the assembly with nanopore reads**

Install Medaka to a new conda environment - https://anaconda.org/bioconda/medaka \
Medaka documentation - https://github.com/nanoporetech/medaka \
Commands used to run Medaka,  \

    conda create env medaka 
    conda create -n medaka -c conda-forge -c bioconda medaka    
    medaka_consensus -i <nanopore reads post filtlong> -d <assembly generated> -o <output directory name> -t <number of threads>
    
**Polishing the assembly with Illumina reads**

Install Pilon - https://github.com/broadinstitute/pilon \
Pilon documentation - https://github.com/broadinstitute/pilon/wiki \
This is an iterative set of commands that would need to be repeated a few times till the number of changes made to the assembly is close to 0. 

    bowtie2-build <input genome> <input genome name>
    bowtie2 -U input_reads.fastq -x <input genome name> --threads 10 -I 0 -X 52 --local --very-sensitive-local | samtools sort > illumina_alignments.bam
    samtools index illumina_alignments.bam
    
    pilon --genome <input genome name> --bam illumina_alignments.bam --output <output genome name> --changes
    rm *.bam *.bam.bai *.bt2
    sed -i 's/_pilon//' <output genome name>
    
### Assembly quality testing

At the end to get the assembly statistics run QUAST in the conda environment genome-assembly geneerated,

    quast.py assembly.fasta -o <assembly-quast-output>

**Test for frameshifts**
Finally long read technologies are error prone, which can annotate lots of pseudogenes due to frameshifts. To determine the quality of the genome assembly on this front, following http://www.opiniomics.org/a-simple-test-for-uncorrected-insertions-and-deletions-indels-in-bacterial-genomes/

    git clone https://github.com/mw55309/ideel.git
    source activate genome-assembly
    conda install -c bioconda prodigal
    conda install -c bioconda diamond
    conda install -c conda-forge r-base
    conda install -c edurand r-cairo
    conda install -c conda-forge libjpeg-turbo
    cd ideel
    #Download uniport database and format UniProt TREMBL and saved as uniprot_trembl.diamond.dmnd

    source activate genome-assembly
    #place the assembly.fasta in the directory ideel/genomes 
    #place uniprot_trembl.diamond.dmnd in ideel/uniprot_trembl.diamond.dmnd
    snakemake --cores 4
    
In the case that there are very few frameshifts in the genome, there will be a tall peak in the graph at x axis 1 suggesting most of the genes in the assemly mapped exactly to proteins to the genes in the database, as shown below. 
<img width="355" alt="image" src="https://user-images.githubusercontent.com/8825721/158039939-a65ea44b-eac3-4831-bb0f-4a2379235244.png">

### Alternative assembly programs to test
In cases where there are only a few fragments as shown above in the graph, its possible these fragments can be resolved by other assemblers. Commands for these other assemblers are provided below

    flye --nano-raw <nanopore reads after QC> --threads 10 --plasmids --out-dir flye-assembly 
    miniasm_and_minipolish.sh <nanopore reads after QC> 10 >miniasm_minipolish_assembly.gfa 
    raven --threads 10 <nanopore reads after QC> >raven_assembly
