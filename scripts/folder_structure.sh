#!/bin/bash

# This script duplicates the folder substructure into a new folder
# The first input is the path to the source folder, and must be in 'name/' format
# The secodn input is the path to the destination folder, which need not yet exist


WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

rsync -a --include '*/' --exclude '*' $1 $2

