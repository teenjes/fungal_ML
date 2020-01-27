#!/bin/bash                                                                                             
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "This script extracts ${RED}$1${NC}"


for file in analysis/Length_Filtered/*/*/length_restricted_reads.fasta; do
	        python scripts/get_align_seqs.py $file $1 -v
	done
