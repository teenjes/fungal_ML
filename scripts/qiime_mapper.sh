#!/bin/bash                                                                                             

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

for file in analysis/Alignment/*/*/1000_reads.fasta; do
		minimap2 -cx map-ont database/sh_refs_qiime_ver8_dynamic_02.02.2019.fasta $file > ${file%1000_reads.fasta}1000_reads_qiime.paf
	done


