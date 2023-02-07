"""
Bacterial genome assembly and annotation

2023, Bhavya Papudeshi 

This is the main snakemake pipeline file.
"""

"""CONFIGFILE
Read in the system default configfile is useful to fill in anything the user has accidentally deleted.
Read in other config files that contain immutable settings that the user is NOT allowed to change.
"""
configfile: os.path.join(workflow.basedir, "config", "config.yaml")

"""ADD FUNCTIONS
If your pipeline needs any Python functions, putting them in a separate file keeps things neat.
"""
#include: "rules/0.functions.smk"

"""PREFLIGHT CHECKS
Validate your inputs, set up directories, parse your config, etc.
"""
include: "rules/1.preflight.smk"

"""TARGETS
Declare your targets, either here, or in a separate file.
"""
include: "rules/2.targets.smk"

"""RULES
Add rules files with the include directive here, or add rules AFTER rule 'all'.
"""
if config['sequencing'] == 'paired':
    include: "rules/qc_qa.smk"
elif config['sequencing'] == 'longread':
    include: "rules/qc_qa_nano.smk"

#add the assembler you choose here, hopefully I have the assembler added as snakemake file already
if config['sequencing'] == 'longread':
    include: "rules/assembly_flye.smk"
elif config['sequencing'] == 'paired':
    include: "rules/assembly_megahit.smk"

#automating picking the phage contigs - circular graphs with high coverage from the assembly results
# Step 1: calculating coverage of all the contigs assembled using coverm 
if config['sequencing'] == 'paired':
    include: "rules/contig_coverage.smk"
    include: "rules/quast.smk"
elif config['sequencing'] == 'longread':
    include: "rules/contig_coverage_nanopore.smk"
    include: "rules/quast_nano.smk"

"""RUN SNAKEMAKE!"""
rule all:
    input:
        allTargets

#rule for just running the assembly, although the default rule all also performs these steps
rule assembly:
    input:
        allTargets
