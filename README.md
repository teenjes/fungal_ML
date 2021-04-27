# Title to be decided

Repository for code involved in filtering DNA sequencing reads, generating machine learning models and linking them together to extract a measure of accuracy at each taxonomic rank, and apply both minimap2 and kraken2 to databases.

## Table of Contents

- database
  - fungal_scpecies_consensus_seqs_44.txt
    - *the gold standard database conaining consensus sequences for each of the 44 species*
- scripts
  - Notebooks
    - kraken2_application.ipynb
      - *jupyter notebook for applying kraken2 to species for accuracy at each taxonomic rank*
    - machine_learning_application.ipynb
      - *jupyter notebook for applying ML to species for accuracy at each taxonomic rank*
    - machine_learning_python_function_gen.ipynb
      - *jupyter notebook for generating python script to creat each ML model*

## Required packages
**Bash**
> minimap2
> kraken2
**Python**
> argparse
> Biopython
> ete3 - NCBITaxa
> json
> keras
> math
> matplotlib
> numpy
> os
> pandas
> random
> sklearn
> subprocess

## 
