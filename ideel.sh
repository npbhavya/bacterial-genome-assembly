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
    #save the input files in the correct directories
    #place the assembly.fasta in the directory ideel/genomes 
    #uniprot_trembl.diamond.dmnd in ideel/uniprot_trembl.diamond.dmnd
    snakemake --cores 4
