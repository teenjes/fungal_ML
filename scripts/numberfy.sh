#!/bin/bash                                                                                                          
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'
PURPLE='\033[1;35m'

echo -e "${GREEN}This script takes each length_restricted fasta file in the ${RED}"$1"${GREEN} directory and converts the base coding to numerical coding for machine learning purposes.

It ${PURPLE} DOES NOT ${GREEN} include padding${WHITE}"


for file in "$1"*/*/length_restricted_reads.fasta; do
	        echo -e "${GREEN} The input file is $file ${WHITE}"
		python ../scripts/numberfy_fasta.py $file
	done
