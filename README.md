# Title to be decided

Repository for code involved in filtering DNA sequencing reads, generating machine learning models and linking them together to extract a measure of accuracy at each taxonomic rank, and apply both minimap2 and kraken2 to databases.

## Table of Contents

- database
  - fungal_scpecies_consensus_seqs_44.txt
    - *the gold standard database containing consensus sequences for each of the 44 species*
- scripts
  - Notebooks
    - kraken2_application.ipynb
      - *jupyter notebook for applying kraken2 to species for accuracy at each taxonomic rank*
    - machine_learning_application.ipynb
      - *jupyter notebook for applying ML to species for accuracy at each taxonomic rank*
    - machine_learning_python_function_gen.ipynb
      - *jupyter notebook for generating python script to creat each ML model*
  - homology_filtering_commands.sh
    - *command for application of minimap2 to original fastq files for homology filtering*
  - kraken2_commands.sh
    - *command for application of kraken2 to original fastq files for homology filtering*
  - length_filtering_and_stats_commands.py
    - *python notebook for applying length filtering and generating summary tables and figures of data*
  - minimap2_commands.sh
    - *command for application of minimap2 to original fastq files for alignment of reads to the gold standard and NCBI databases*

## Required packages
**Bash**
> minimap2 <br>
> kraken2 <br>
<br>

**Python**
> argparse <br>
> Biopython <br>
> ete3 - NCBITaxa <br>
> json <br>
> keras <br>
> math <br>
> matplotlib <br>
> numpy <br>
> os <br>
> pandas <br>
> random <br>
> sklearn <br>
> subprocess <br>

## Credits

Tavish Eenjes <br>
Yiheng Hu
