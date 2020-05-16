#!/bin/bash

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

rm database/custom_database_100.fasta
for file in analysis/Consensus/*/*/clean_consensus_100.fasta; do
	        python scripts/make_database_100.py $file -v
	done
