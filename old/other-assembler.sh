#!/bin/bash

#SBATCH --job-name=assemblers
#SBATCH --mail-user=<email>
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out.txt
#SBATCH --error=%x-%j.err.txt
#SBATCH --time=1-0
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G

module load miniconda/3.0 
source activate genome

cd <filepath>

flye --nano-raw <nanopore reads after QC> --threads 10 --plasmids --out-dir flye-assembly 
miniasm_and_minipolish.sh <nanopore reads after QC> 10 >miniasm_minipolish_assembly.gfa 
raven --threads 10 <nanopore reads after QC> >raven_assembly
