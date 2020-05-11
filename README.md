# 2019/2020 Honours Project
## Developing Tools for Fungal Pathogen Identification
### Tavish Eenjes


**This repository contains all data and code required for the project**

Contents
========
- raw_data (Directory)
	- *Contains the README files for the primary dataset sequenced by Yiheng*
		- error_profile_README.txt - primary dataset
	
- analysis (Directory)
	- *Contains all directories and files related to the data pre-processing, filtering and analysis*
		- Python_Processing (Directory)
			- homology (Directory)
				- homology_genus_info.txt
					- *Text file detailing the sources for all sequences used for the generation of the homology_genus.fasta database*
				- homology_genus.fasta
					- *Database file containing the sequences used for frDNA homology filtering*
- scripts (Directory)
	- deepbinner_script.sh
		- *Script for taking an input directory containing mixed fast5 reads and processing them through Deepbinner, specifying an output directory for demultiplexed single-read fast5 files*
	- guppy_script.sh
                - *Script for taking all folders within an input directory containing demultiplexed multi-read fast5 files and basecalling them using the Guppy basecaller, specifying a output directory for basecalled fastq files*
	- homology_map.sh
		- *Script for taking all folders within an input directory and applying the minimap2 program to map reads against a custom-curated database formed of 29 sequences (analysis/Python_Processing/homology/homology_genus_info.txt) for assessing frDNA homology*
	- ontfast5api.sh
		- *Script for taking an input directory containing single-read fast5 files and processing them through the ont_fast5_api software from Oxford Nanopore, specifying an output directory for multi-read fast5 files*
	- tensorflow_test.py
		- *Script for testing whether TensorFlow, required for Deepbinner, is being run on CPU or GPU*

