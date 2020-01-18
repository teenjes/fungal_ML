#!/bin/bash                                                                                             
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


echo -e "This script takes the consensus_100.fasta file from each species and runs the EMBOSS cons software on it to generate a consensus sequence, then subsequently runs a cleanup python script to remove any 'N's or 'n's from the resultant consensus"

for file in analysis/Consensus/*/*/consensus_100.fasta; do
		cons -sequence $file -outseq ${file%consensus_100.fasta}/cons_consensus.fasta -identity 50 -name ${file%consensus_100.fasta}
		python scripts/cleanup.py ${file%consensus_100.fasta}/cons_consensus.fasta -v
done



#for i in $(seq 01 100); do
#	        cons -sequence analysis/Consensus/20171103_FAH15473/barcode02/consensus_100.fasta -outseq analysis/Consensus/20171103_FAH15473/barcode02/test$i.fasta -identity $i -name test$i
		        #python ???.py ???folder???
#		done
