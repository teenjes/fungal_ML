#!/bin/bash

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

echo -e "${GREEN}This script takes all files within the input directory ${RED}"$1"${GREEN} and converts them to multi_fast5 files, depositing them in the output path, ${RED}"$2"${WHITE}"

for d in "$1"/*/; do
	name="$(basename -- $d)"
	input="$d"
	output="$2/$name"
	echo "input is $input"
	echo "output is $output"
	echo -e "${RED} BEGIN ${WHITE}"
	single_to_multi_fast5 -i $input -s $output
	echo -e "${RED} COMPLETE ${WHITE}"
done
