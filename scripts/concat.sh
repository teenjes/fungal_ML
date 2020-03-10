#!/bin/bash                                                                                                            






WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "${GREEN}This script takes all directories within the input directory ${RED}"$1"${GREEN} and finds the barcode directories containing fastq files. It then runs a cat command to merge these files and save them to an output directory ${RED}"$2"${WHITE}"

for d in "$1"/*/; do
	name="$(basename -- $d)"
        input="$d"
	output="$2/$name"
	echo "input is $input"
	echo "output is $output"
	if [ ! -d "$output" ]; then
		mkdir "$output"
		echo -e "${RED}This directory was created${WHITE}"
	fi
	cat $input*.fastq > $output/merged.fastq
done
