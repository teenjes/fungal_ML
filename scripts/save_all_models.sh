#!/bin/bash                                                                                             
WHITE='\033[1;37m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

while read tax_name name; do
	echo "tax_name is: $tax_name"
	echo "name is: $name"
	python scripts/machine_learning.py analysis/Stats/reference_dataframe.csv analysis/ -r $tax_name -n $name -c 15000 -v
done < $1
