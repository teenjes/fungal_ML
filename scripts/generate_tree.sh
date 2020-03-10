#!/bin/bash
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

# 100
for file in analysis/Length_Restricted/*/*/length_restricted_reads.fasta; do
	                python scripts/get_cons_ids3.py $file -v
			                        rm *.logfile
					done


for file in analysis/Consensus/*/*/for_consensus_100.fasta; do
	                echo -e "${GREEN} The input file is $file ${WHITE}"
			                        STR=$file
						                                starter=${file%/*}
										                                        ender='for_consensus_100.aln'
															                                                comb=$starter/$ender
																					                                                        muscle -in $file -out $comb -maxiters 8
																												                                                                echo -e "${RED} The output file is $comb ${WHITE}"
																																				                                                        done
																																											for file in analysis/Consensus/*/*/for_consensus_100.aln; do
																																												                        cons -sequence $file -outseq ${file%for_consensus_100.aln}/cons_consensus_101.fasta -name ${file%for_consensus_100.fasta}
																																															                                        python scripts/cleanup_101.py ${file%for_consensus_100.aln}/cons_consensus_101.fasta -v
																																																				                                done
																																																								rm database/custom_database_101.fasta
																																																								for file in analysis/Consensus/*/*/clean_consensus_101.fasta; do
																																																									                python scripts/make_database_101.py $file -v
																																																										done
																																																										python scripts/phylogenetic_naming_101.py database/custom_database_101.fasta -v

