#!/bin/bash                                                                                                          
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "${GREEN}This script works through all files (using those in the ${RED}"$1"${GREEN} directory for reference) and runs the summary_statistics python script to get summary statistics and figures for each file ${WHITE}"

for i in $(seq 01 100); do 
	cons -sequence analysis/Consensus/20171103_FAH15473/barcode02/consensus_100.fasta -outseq analysis/Consensus/20171103_FAH15473/barcode02/test$i.fasta -identity $i -name test$i
	#python ???.py ???folder???
done
