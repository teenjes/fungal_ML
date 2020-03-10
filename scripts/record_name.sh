#!/bin/bash




echo "Run from the analysis/Fungi folder"


WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


for d in */*/*/*/*/*/r*; do
	#echo "$d"
	python ../../scripts/record_name.py "$d" -v
done

for d in */*/*/*/*/*/*/r*; do
	#echo "$d"
	python ../../scripts/record_name.py "$d" -v
done
