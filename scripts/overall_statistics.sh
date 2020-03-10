#!/bin/bash
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

rm analysis/Stats/read_number.csv
for file in analysis/Concatenated/*/*/merged.fastq; do
	        python scripts/overall_statistics.py $file -v
	done
