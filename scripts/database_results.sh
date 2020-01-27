#!/bin/bash                                                                                             
WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

for file in analysis/Alignment/*/*/minimap_result.paf; do
	        python scripts/1000_minimap_result.py $file -v
	done
