#!/bin/bash

# This script is designed to be run such that deepbinner is performed on a dataset to demultiplex reads

in_loc="$1"
out_loc="$2"

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

echo -e "${GREEN}This script takes all folders within the input directory ${RED}"$1"${GREEN} and runs Deepbinner of these files to demultiplex the fast5 files within these folders and output them into barcode directories within the output directory ${RED}"$2"${WHITE}"

"
for d in "$1"*/ ; do
	deepbinner realtime --in_dir "$d" --out_dir "$2" --native
done"
