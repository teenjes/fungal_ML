#!/bin/bash
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

for file in analysis/Consensus/*/*/; do
	cat $file/clean_consensus.fasta <(echo) $file/clean_consensus_100.fasta <(echo) $file/clean_consensus_101.fasta <(echo) $file/clean_consensus.fasta <(echo) $file/clean_consensus_100.fasta <(echo) $file/clean_consensus_101.fasta <(echo) $file/clean_consensus.fasta <(echo) $file/clean_consensus_100.fasta <(echo) $file/clean_consensus_101.fasta <(echo) $file/clean_consensus.fasta <(echo) $file/clean_consensus_100.fasta <(echo) $file/clean_consensus_101.fasta > $file/combined_clean_consensus.fasta
done
