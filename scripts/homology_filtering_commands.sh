#!/bin/bash


echo -e "This script works through all folders within the input directory "$1" and uses the minimap2 program to map reads against a custom-curated database formed of 29 sequences for assessing frDNA homology"

# for files in the original unfiltered fastq, apply minimap2 for homology filtering
for d in "$1"/Concatenated; do
	for f in "$d"/*; do
		echo "$f"
		remove='analysis/Concatenated'
		output=${f//$remove}
		echo "$output"
		minimap2 -x map-ont analysis/Python_Processing/homology/homology_genus.fasta "$f"/merged.fastq > analysis/Python_Processing"$output"/combined_test.paf
done
done
