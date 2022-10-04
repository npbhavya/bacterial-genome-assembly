import os
from setuptools import setup

def get_version():
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'bacterial_genome', 'bacterial_genome.VERSION')) as f:
        return f.readline().strip()

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
 name='bacterial-genome',
 description="Assembling bacterial genomes and annotations ",
 version=get_version(),
 author="Bhavya Papudeshi",
 author_email="npbhavya13@gmail.com",
 py_modules=['bacterial_genome'],
 install_requires=["snakemake==7.14.0",
                   "pyyaml==6.0",
                   "Click==8.1.3"],
 entry_points={
  'console_scripts': [
    'bacterial_genome=bacterial_genome.__main__:main'
  ]},
 include_package_data=True,
)
