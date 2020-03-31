#!/bin/bash
# This should be run in the ANALYSIS folder

WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

for file in analysis/database_mapping/custom/*.paf; do
	python scripts/custom_minimap_result.py $file
done
