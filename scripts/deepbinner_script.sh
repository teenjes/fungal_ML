#!/bin/bash

# This script is designed to be run such that deepbinner is performed on a dataset to demultiplex reads
# It requires an input directory that contains directories with all reads from a single MinION flowcell - these subdirectories are examined one-by-one by the script
# It requires an output directory as a location for the deposit of demultiplexed reads. This should be labelled with the flowcell run

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

echo -e "${GREEN}This script takes all folders within the input directory ${RED}"$1"${GREEN} and runs Deepbinner of these files to demultiplex the fast5 files within these folders and output them into barcode directories within the output directory ${RED}"$2"${WHITE}"


for d in "$1"/*/; do
	echo "$d"
	echo "$2"
	deepbinner realtime --in_dir "$d" --out_dir "$2" --native
done

