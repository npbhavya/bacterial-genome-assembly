"""

Rules for quality control and quality assurance - Illumina paired end reads 
"""

rule quast_nano:
    input:
        contigs =os.path.join(ASSEMBLY, "{sample}-flye", "assembly.fasta"),
    output:
        stats = os.path.join(ASSEMBLY, "{sample}_flye_quast/report.txt"),
    conda: "../envs/quast.yaml"
    params:
        o = os.path.join(ASSEMBLY, "{sample}_flye_quast")
    log:
        os.path.join(logs, "quast_{sample}.log")
    shell:
        """
        quast.py {input.contigs} -o {params.o} 2> {log}
        """
