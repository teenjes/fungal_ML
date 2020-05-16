#!/bin/bash                                                                                                         
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


for file in "$1"*/*/for_consensus.fasta; do
	        echo -e "${GREEN} The input file is $file ${WHITE}"
		        STR=$file
			        starter=${file%/*}
				        ender='for_consensus.aln'
					        comb=$starter/$ender
						        muscle -in $file -out $comb -maxiters 8
							        echo -e "${RED} The output file is $comb ${WHITE}"
							done
