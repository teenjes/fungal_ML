#!/bin/bash







WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "${GREEN}This script takes all folders within the input directory ${RED}"$1"${GREEN} and runs the Guppy basecaller r9.4.1 on these files to decode the base-pair data and output them into barcode directories within the output directory ${RED}"$2"${WHITE}"

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
	guppy_basecaller -i "$input" -s "$output" -c "dna_r9.4.1_450bps_hac.cfg" --device "auto"
done
