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

cd genomes

export LC_ALL=en_US.UTF-8
filtlong --min_length 1000 --keep_percent 95 <fastq files from guppy> >> <genome name>_filtlong.fastq
prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 1 -derep 1 -out_format 0 -trim_tail_left 5 -trim_tail_right 5 -ns_max_n 5  -trim_qual_type min -trim_qual_left 30 \
			-trim_qual_right 30 -trim_qual_window 10 -threads 10 -out_good <genome name>_prinseq_good1 -out_single <genome name>_prinseq_single1 -out_bad <genome name>_prinseq_bad1 \
			-out_good2 <genome name>_prinseq_good2 -out_single2 <genome name>_prinseq_single2 -out_bad2 <genome name>_prinseq_bad2 \
			-fastq <illunina forward file> -fastq2 <illumina reverse file>
      
unicycler -1 <genome name>_prinseq_good -2 <genome name>_prinseq_good2 -l <genome name>_filtlong.fastq -o unicycler_assembly --no_correct"
