#!/bin/bash                                                                                             

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

for file in analysis/Alignment/*/*/1000_reads.fasta; do
	minimap2 -cx map-ont $1 $file > ${file%1000_reads.fasta}minimap_result.paf
done

