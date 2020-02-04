#!/bin/bash
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "This script takes the consensus2_1000.fasta file from each species and runs the EMBOSS cons software on it to generate a consensus sequence, then subsequently runs a cleanup python script to remove any 'N's or 'n's from the resultant consensus"

for file in analysis/Consensus/20171103_FAH15473/*/consensus2_1000.fasta; do
	                cons -sequence $file -outseq ${file%consensus2_1000.fasta}/cons_consensus_1000.fasta -name ${file%consensus2_1000.fasta}
			                python scripts/cleanup.py ${file%consensus2_1000.fasta}/cons_consensus_1000.fasta -v
				done
