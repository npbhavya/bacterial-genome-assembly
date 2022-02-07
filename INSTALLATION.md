# Dependecies

- miniconda https://docs.conda.io/en/latest/miniconda.html
- Jupyter notebook https://jupyter.org/
- Bandage plots https://rrwick.github.io/Bandage/
- unicycler https://anaconda.org/bioconda/unicycler
- prinseq++ https://anaconda.org/bioconda/prinseq-plus-plus
- filtlong https://anaconda.org/bioconda/filtlong
- quast https://anaconda.org/bioconda/quast
- snakemake https://anaconda.org/bioconda/snakemake


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
  `pip install quast`
