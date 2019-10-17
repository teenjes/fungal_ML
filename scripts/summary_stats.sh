#!/bin/bash

# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "${GREEN}This script works through all files (using those in the ${RED}"$1"${GREEN} directory for reference) and runs the summary_statistics python script to get summary statistics and figures for each file ${WHITE}"


for subdirectory in "$1"*; do
	for folder in "$subdirectory"/*; do
		for file in "$folder"/*; do
			remove1='Concatenated/'
			remove2='merged.fastq'
			output1=${file//$remove1/}
			output2=${output1//$remove2/}
			python ../scripts/summary_statistics.py "$file" Python_Processing/"$output2" Stats/"$output2" -v
		done
	done
done

