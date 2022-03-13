#!/bin/bash
            #SBATCH --job-name=pilon
            #SBATCH --mail-user=email
            #SBATCH --mail-type=ALL
            #SBATCH --output=%x-%j.out.txt
            #SBATCH --error=%x-%j.err.txt
            #SBATCH --time=0-4
            #SBATCH --ntasks=10
            #SBATCH --cpus-per-task=1
            #SBATCH --mem=200G

            cd <file path>

            module add bowtie/2.3.5.1
            module add samtools

            bowtie2-build <input genome> <input genome name>
            bowtie2 -U input_reads.fastq -x <input genome name> --threads 10 -I 0 -X 52 --local --very-sensitive-local | samtools sort > illumina_alignments.bam
            samtools index illumina_alignments.bam

            module add miniconda/3.0
            source activate genome
            pilon --genome <input genome name> --bam illumina_alignments.bam --output <output genome name> --changes
            #rm *.bam *.bam.bai *.bt2
            sed -i 's/_pilon//' <output genome name>
