#!/bin/bash                                                                                             # This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "${GREEN}This script works through all files (using those in the ${RED}"$1"${GREEN} directory for reference) and runs the summary_statistics python script to get summary statistics and figures for each file ${WHITE}"

# 100
for file in "$1"*/*/sequencing_summary.txt; do
	                python ../scripts/for_consensus_high_qscore_100.py $file -v
			                        rm *.logfile
						                done
