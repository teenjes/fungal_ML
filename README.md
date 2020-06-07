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
	- cleanup_100.py
		- *Scrpt for removing non-nucleotide characters from a consensus sequence*
	- create_consensus_100.sh
		- *Script for applying Muscle to generate an alignment file for use in consensus generation*
	- deepbinner_script.sh
		- *Script for taking an input directory containing mixed fast5 reads and processing them through Deepbinner, specifying an output directory for demultiplexed single-read fast5 files*
	- emboss_cons_100.sh
		- *Script for creating a consensus sequence from an alignment file, then subsequently applying the cleanup_100.py script to remove non-nucleotide characters eg. 'n' or 'N'*
	- get_cons_ids.py
		- *Script that applies Mothur to remove primer sequences and subsamples 100 reads for use in consensus generation*
	- guppy_script.sh
                - *Script for taking all folders within an input directory containing demultiplexed multi-read fast5 files and basecalling them using the Guppy basecaller, specifying a output directory for basecalled fastq files*
	- homology_map.sh
		- *Script for taking all folders within an input directory and applying the minimap2 program to map reads against a custom-curated database formed of 29 sequences (analysis/Python_Processing/homology/homology_genus_info.txt) for assessing frDNA homology*
	- make_database_100.sh
		- *Script to collectall consensus sequences in the one file*
	- ontfast5api.sh
		- *Script for taking an input directory containing single-read fast5 files and processing them through the ont_fast5_api software from Oxford Nanopore, specifying an output directory for multi-read fast5 files*
	- summary_statistics.py
		- *Python script that utilises the restricted list of read IDs generated by homology filtering, and applied further length filtering to output the length-filtered output fasta file for use in further analysis. This script also then generates and saves figures showing the before and after filtering spread of reads*
	- summary_stats.sh
		- *Bash script that loops an automates the summary_statistics python script across each sample present in the input directory*
	- tensorflow_test.py
		- *Script for testing whether TensorFlow, required for Deepbinner, is being run on CPU or GPU*
	- Notebooks (Directory)
		- Alignment-Based_Workbook.ipynb
			- *Python notebook for using minimap2 to map reads against the Consensus and Qiime databases*
		- Application_alignment.ipynb
			- *Python notebook for subsampling the wheat and synthetic mock community datasets, then applying the alignment-based technique to both datasets for each of the Consensus and Qiime databases*
		- Crossmapping.ipynb
			- *Python notebook for assessing the relationship between the coss-mapping fraction for reads against the Consensus database and the genetic distance between species*
		- machine_learning_script.ipynb
			- *Python notebook for the generation of machine learning models, and for the use of those models when applied to the synthetic mock community and wheat datasets*

