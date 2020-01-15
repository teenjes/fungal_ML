#!/bin/bash                                                                                                          # This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "${GREEN}This script removes any existing scores.txt file. and then accepts test*.fasta files in the ${RED}"$1"${GREEN} directory, as consensus sequences to compare to 100 sequences randomly subsampled from the full population ${WHITE}"

rm $1scores.txt

for file in $1test*.fasta; do
	python scripts/consensus_test.py $file -v
done
