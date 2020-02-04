#!/bin/bash                                                                                                            
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "${GREEN}This script works through all files (using those in the ${RED}"$1"${GREEN} directory for reference) and runs the summary_statistics python script to get summary statistics and figures for each file ${WHITE}"


for file in "$1"*/*/for_consensus_1000.fasta; do
	echo -e "${GREEN} The input file is $file ${WHITE}"
	STR=$file
	starter=${file%/*}
        ender='for_consensus_1000.aln'
	comb=$starter/$ender
	muscle -in $file -out $comb -maxiters 2
	echo -e "${RED} The output file is $comb ${WHITE}"
done
