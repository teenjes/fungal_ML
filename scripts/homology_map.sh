#!/bin/bash


WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

echo -e "${GREEN}This script works through all folders within the input directory ${RED}"$1"${GREEN} and uses the minimap2 program to map reads against a custom-curated database formed of 29 sequences for assessing frDNA homology ${WHITE}"


for d in "$1"/Concatenated; do
	for f in "$d"/*; do
		echo "$f"
		remove='analysis/Concatenated'
		output=${f//$remove}
		echo "$output"
		minimap2 -x map-ont analysis/Python_Processing/homology/homology_genus.fasta "$f"/merged.fastq > analysis/Python_Processing"$output"/combined_test.paf
done
done
