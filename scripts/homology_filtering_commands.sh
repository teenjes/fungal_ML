#!/bin/bash

### Requires minimap2 and ONT reads


echo -e "This script takes an input fastq file "$1" and uses the minimap2 program to map reads against a custom-curated database ($3), saving to output directory "$2""

# for files in the original unfiltered fastq, apply minimap2 for homology filtering

# input = $1
# output = $2
# custom database is at $3



minimap2 -x map-ont $3 $1 > "$2"/mapped.paf
