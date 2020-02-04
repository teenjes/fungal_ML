#!/bin/bash
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "This script takes the consensus2_1000.fasta file from each species and runs the EMBOSS cons software on it to generate a consensus sequence, then subsequently runs a cleanup python script to remove any 'N's or 'n's from the resultant consensus"

for file in analysis/Consensus/*/*/for_consensus_100.aln; do
	                        cons -sequence $file -outseq ${file%for_consensus_100.aln}/cons_consensus_101.fasta -name ${file%for_consensus_100.fasta}
				                                        python scripts/cleanup_101.py ${file%for_consensus_100.aln}/cons_consensus_101.fasta -v
									                                done
