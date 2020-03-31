#!/bin/bash                                                                                             
WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

for file in analysis/database_mapping/custom/*.fasta; do
	        minimap2 -cx map-ont $1 $file > ${file%.fasta}.paf
	done
