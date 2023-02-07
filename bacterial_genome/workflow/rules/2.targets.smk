"""
Declare your targets here!
A separate file is ideal if you have lots of target files to create, or need some python logic to determine
the targets to declare. This example shows targets that are dependent on the input file type.
"""
allTargets =[]


#QC_QA
if config['sequencing'] == 'paired':
    allTargets.append(expand(os.path.join(QCDIR, "{sample}_good_out_R1.fastq"), sample=SAMPLES))
    allTargets.append(expand(os.path.join(QCDIR, "{sample}_good_out_R2.fastq"), sample=SAMPLES))
elif config['sequencing'] == 'longread':
    allTargets.append(expand(os.path.join(QCDIR, "{sample}-filtlong.fastq"), sample=SAMPLES))

#assembly
if config['sequencing'] == 'longread':
    allTargets.append(expand(os.path.join(ASSEMBLY, "{sample}-flye", "assembly.fasta"), sample=SAMPLES))
    allTargets.append(expand(os.path.join(ASSEMBLY, "{sample}-flye", "assembly_graph.gfa"), sample=SAMPLES))
    allTargets.append(expand(os.path.join(ASSEMBLY, "{sample}-flye", "assembly_info.txt"), sample=SAMPLES))
elif config['sequencing'] == 'paired':
    allTargets.append(expand(os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.contigs.fa"), sample=SAMPLES))
    allTargets.append(expand(os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.fastg"), sample=SAMPLES))

#assembly stats
if config['sequencing'] == 'paired':
    allTargets.append(expand(os.path.join(ASSEMBLY, "{sample}_megahit_quast/report.txt"), sample=SAMPLES))
    allTargets.append(expand(os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}-contigs.tsv"), sample=SAMPLES))
elif config['sequencing'] == 'longread':
    allTargets.append(expand(os.path.join(ASSEMBLY, "{sample}_flye_quast/report.txt"), sample=SAMPLES))
    allTargets.append(expand(os.path.join(ASSEMBLY, "{sample}-flye", "{sample}-contigs.tsv"), sample=SAMPLES))

