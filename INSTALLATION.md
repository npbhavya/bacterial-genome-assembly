# Dependecies

- miniconda https://docs.conda.io/en/latest/miniconda.html
- Jupyter notebook https://jupyter.org/
- Bandage plots https://rrwick.github.io/Bandage/
- prinseq++ https://anaconda.org/bioconda/prinseq-plus-plus
- filtlong https://anaconda.org/bioconda/filtlong
- snakemake https://anaconda.org/bioconda/snakemake
- unicycler https://anaconda.org/bioconda/unicycler
- flye https://anaconda.org/bioconda/flye
- raven https://anaconda.org/conda-forge/raven
- miniasm https://anaconda.org/bioconda/miniasm
- quast https://anaconda.org/bioconda/quast

# INSTALL 

Create a new conda enviornment 

  `conda env create -n genome-assembly`

Activate the environment 

  `source activate genome-assembly`
  
Install dependencies 

  `conda install -c bioconda prinseq` \
  `conda install -c bioconda unicycler` \
  `conda install -c bioconda flye` \
  `conda install -c bioconda filtlong` \
  `conda install -c conda-forge mamba` \
  `conda install -c bioconda snakemake` \
  `conda install -c conda-forge raven` \
  `conda install -c bioconda miniasm` \
  `pip install quast`
  
  Visualization of the assemblies 
  - Install banadage https://rrwick.github.io/Bandage/ on laptop
