#!/bin/bash
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

for file in analysis/Length_Filtered/*/*/length_restricted_reads.fasta; do
	cp $file ${file//length_restricted_reads.fasta}length_restricted_for_use.fasta
done       
