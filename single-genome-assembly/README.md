# Running genome assembly on deepthought Flinders University cluster
## For a single genome file sequenced on Nanopore MinION platform

Documentation on deepthought - https://deepthoughtdocs.flinders.edu.au/en/latest/index.html

## fast5 basecalling using guppy-basecaller 

**Download** 

Available on Nanopore community website. Download the CPU or GPU version

        tar -xvzf <binary file>

**Run as script**

Save the below script as a job file, changing email, and filepaths. 

        #!/bin/bash

        #SBATCH --job-name=guppy_basecaller
        #SBATCH --mail-user=email
        #SBATCH --mail-type=ALL
        #SBATCH --output=%x-%j.out.txt
        #SBATCH --error=%x-%j.err.txt
        #SBATCH --time=1-0
        #SBATCH --ntasks=1
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=50G
        #SBATCH --gres=gpu:1
        #SBATCH --partition=hpc_gpu

        module load cuda10.0/nsight/10.0.130
        export PATH=~/ont-guppy/bin:$PATH

        cd <path to run the script and find the files>
        guppy_basecaller --input_path <directory with fast5 files> --save_path <output path> --flowcell FLO-             FLG001 --kit SQK-RBK004 --device auto

To submit the script 

        sbatch <job script.sh>
        
**Output files**

fastq files saved in output path defined with multiple fastq files.

        cat <output path>/*.fastq >> sample_basecalled.fastq 
        
This is the input fastq files for genome assembly

## Genome assembly steps

**loading the environment**

    module load Miniconda/3.0

**create a new conda environment**

    conda create -y -n genome-assembly
    source activate genome-assembly
    conda install -c bioconda unicycler
    conda install -c bioconda filtlong
    conda install -c conda-forge mamba
    conda install -c conda-forge -c bioconda
  
**Input files** 

Fastq read file (post basecalling)

**Steps run in the job script below**
- QC using filtlong - https://github.com/rrwick/Filtlong
- Assembly using Unicycler - https://github.com/rrwick/Unicycler

**wrote a job script and submitted to SLURM job scedule** 

    #!/bin/bash
    #SBATCH --job-name=assembly
    #SBATCH --mail-user=email
    #SBATCH --mail-type=ALL
    #SBATCH --output=%x-%j.out.txt
    #SBATCH --error=%x-%j.err.txt
    #SBATCH --time=0-24
    #SBATCH --ntasks=10
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=250G

    module load miniconda/3.0 
    source activate genome-assembly

    cd Acinetobacter

    export LC_ALL=en_US.UTF-8
    filtlong --min_length 1000 --keep_percent 95 recalled.fastq >Acinetobacter.fastq 
    unicycler -l Acinetobacter.fastq -o assembly --no_correct 
    
**Output** 

The output if the script ran successfully, will generate an output directory with assembly.fasta and assembly.gfa files.

## Visualization using Bandage plot
Download the assembly.gfa file- output from Unicycler output locally to visualize on Bandage- https://rrwick.github.io/Bandage/
Upload the assembly.gfa file to visulaize the genome assembled


## Assembly quality
Ideally, the assembly will generate complet, circular chromosomes along with complete, circular plasmid sequences if present in the genome sequenced. 
Number of contigs =Number of chromsomes in the genome +number of plasmid present in the genome 
Total length=Total genome length of species sequenced 

Finally long read technologies are error prone, which can annotate lots of pseudogenes due to frameshifts. To determine the quality of the genome assembly on this front, following http://www.opiniomics.org/a-simple-test-for-uncorrected-insertions-and-deletions-indels-in-bacterial-genomes/

**Setting up the script to test frameshifts

    git clone https://github.com/mw55309/ideel.git
    source activate genome-assembly
    conda install -c bioconda prodigal
    conda install -c bioconda diamond
    conda install -c conda-forge r-base
    conda install -c edurand r-cairo
    conda install -c conda-forge libjpeg-turbo
    cd ideel
    #Download uniport database and format UniProt TREMBL and saved as uniprot_trembl.diamond.dmnd

**Running the scripts**

**Input files** - place the assembly.fasta in the directory ideel/genomes
**Database** - place uniprot_trembl.diamond.dmnd in ideel/uniprot_trembl.diamond.dmnd

**Setting up job script to submit to SLURM job scheduler**

        #!/bin/bash
        #SBATCH --job-name=frameshift_test
        #SBATCH --mail-user=email
        #SBATCH --mail-type=ALL
        #SBATCH --output=ideel/%x-%j.out.txt
        #SBATCH --error=ideel/%x-%j.err.txt
        #SBATCH --time=4:00:00
        #SBATCH --ntasks=6
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=50G

        module load miniconda/3.0 
        source activate genome-assembly

        cd ideel
        snakemake --cores 4
    
**Output**

In ideel/hists/<output filename>.png plots

Ideally - Should have a lot of proteins with 1:1 ratio on the histograms (one tall peek)and no other peaks at all
The other peaks represent proteins that were only partially aligned against the proteins in the reference database, and can be pseudogenes or genes with frameshifts. Detailed explanation in the blog mentioned above.





