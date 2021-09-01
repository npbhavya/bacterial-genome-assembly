# Snakemake workflow

## Input file
Input files are nanopore reads in fast5 format, save the fast5 files to the fast5 directory

## Conda env
The workflow should be run in conda environment "genome-assembly"

	source activate genome-assembly

Next run the snakemake workflow, the workflow starts with 
- fast5 to fastq conversion
- quality control to remove short reads and low quality reads 
- assembly using unicycler, flye, raven and miniasm 
- viralverify to identify viral reads from each of the assemblies
- saving the viral reads to the "viral_contigs" directory for each of the assembly

	snakemake -s snakefile-1 --cores 4 

The cores is set to 4 since guppy is run on gpu nodes which runs only if there are 4 samples running at once. 
After this step though, the rest of the steps are run on cpu, so they can be run with higher number of cores.

## Manual steps

Next, fasta files in "viral_contigs" that are empty are removed

	for f in viral_contigs/* ; do  if [ -s "$f" ]; then  echo "not empty"; else rm -rf "$f" ; fi; done 


Manualy running trycycler, https://github.com/rrwick/Trycycler/wiki/Clustering-contigs onwards till generating consensus assembly
Save the consensus assembly in viral_contigs directory at the end

## Continue snakemake workflow
Run the next snakefile to calculate checkv completeness 
	
	snakemake -s snakefile-2 --cores 10

## Manual steps
Take a look at the output file, quality_summary_all.tsv to pick the assembler that prcovided
- if trycycler, pick this one 
- the highest completeness 
- most closest to the genome length to predicted genome length or the longest
- look at genome quality, number of genes and viral genes


Saved the final assembly - with representation of one assembly per genome 

### Polishing the assemblies 

Medaka polishing, done manually
	
	source activate medaka
	medaka_consensus -i filtlong/Lysate_Run_Bb18-filtlong.fastq -d select_one_assembler/Lysate_Run_Bb18-flye.fasta -o medaka/Lysate_Run_Bb18-medaka -m r941_min_sup_g507  -t 10 


Pilon polishing if Illumina reads available, Bu20, Bu07, Bu18, Bu19. Followd the pilon tutorial in https://github.com/rrwick/Trycycler/wiki/Polishing-after-Trycycler
Saved the output to directory "polishing"/*.fasta

## rerunning checkv again on the polished genomes 
Running snakemake, checkv output

	source activate genome-assembly
	snakemake -s snakemake-3 --cores 10


This runs checkv on polished genomes
rename the contigs
split the assemblies if they have multifasta
runs fasANI on the final set
runs mash on the final set


## Final genomes, manual steps
Removing replicate genomes if ANI are really high (>99%) and mash distances are really low (<0.005)


## Annotation steps

	snakemake -s snakemake-annot --cores 10

- First start with making sure all the genomes start with terminase gene,
- reverse complemented to have all genomes on the smae orientation
- run phanotate or genotate

