#!/bin/bash                                                                                             
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'




for i in $(seq 0 99); do
	minimap2 analysis/Consensus/20171103_FAH15473/barcode02/test$i_new.fasta analysis/Alignment/20171103_FAH15473/barcode02/100_reads.fasta > analysis/Consensus/20171103_FAH15473/barcode02/test$i.paf        
	
		done
