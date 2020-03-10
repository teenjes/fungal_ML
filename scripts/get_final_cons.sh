#!/bin/bash

# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "${GREEN}This script works through all files (using those in the ${RED}"$1"${GREEN} directory for reference) and runs the summary_statistics python script to get summary statistics and figures for each file ${WHITE}"

# 1000
for file in "$1"*/*/length_restricted_reads.fasta; do
	                python ../scripts/get_final_cons.py $file -v
			                        rm *.logfile
						                done
