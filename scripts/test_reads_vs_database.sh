#!/bin/bash
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "This script takes the consensus2_1000.fasta file from each species and runs the EMBOSS cons software on it to generate a consensus sequence, then subsequently runs a cleanup python script to remove any 'N's or 'n's from the resultant consensus"


for file in analysis/Alignment/*/*/1000_reads.fasta; do
	                        minimap2 database/sh_refs_qiime_ver8_dynamic_02.02.2019.fasta $file > analysis/Alignment/${file:19:27}/1000_reads_map.paf
				                                python scripts/test_reads_vs_database.py analysis/Alignment/${file:19:27}/1000_reads_map.paf
								                        done
