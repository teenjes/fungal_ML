
"""
The goal of this program is to examine the distribution of reads within
each file, and for the files generated after homology analysis.
The program will generate summary statistics for the result of the homology
analysis and save figures illustrating the read distribution for the
frDNA reads
"""

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord as SR
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import csv


WHITE='\033[1;37m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'
PURPLE='\033[1;35m'

parser = argparse.ArgumentParser(description="""
The goal of this program is to examine the distribution of reads within
each file, and for the files generated after homology analysis.
The program will generate summary statistics for the result of the homology
analysis and save figures illustrating the read distribution for the
frDNA reads
""")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
parser.add_argument("full_file")
args = parser.parse_args()
print('\033[0;35m'+'START'+'\033[1;37m')
full_file = args.full_file
print(full_file)
homology_folder = args.full_file[:9]+'Python_Processing'+args.full_file[21:-12]+'combined_test.paf'
print(homology_folder)
length_folder = args.full_file[:9]+'Length_Filtered'+args.full_file[21:-12]+'length_restricted_reads.fasta'
print(length_folder)

# Load the full file containing all reads for this barcode
full_file_dict = SeqIO.to_dict(SeqIO.parse(full_file, "fastq"))

if args.verbose:
    print('\033[0;34m' + "Loaded " + full_file + '\033[1;37m')

    
# Open the original merged.fastq file
# Extract the number of reads contained within this file
full_lengths = []
for key in full_file_dict:
    full_lengths.append(len(full_file_dict[key].seq))
full_lengths_len = len(full_file_dict)
print(full_lengths_len)


# Import the PAF file resulting from the minimap2 homology filtering
homology_paf = pd.read_csv(homology_folder, sep='\t', header=None, engine='python')

if args.verbose:
    print('\033[0;34m' + "Loaded " + homology_folder + '\033[1;37m')
          
# Create new dictionary using only keys present in paf file
# Extract the number of reads contained within this file
homology_dict = {}
for key in homology_paf[0].unique():
    homology_dict[key] = full_file_dict[key]

homology_lengths = []
for key in homology_dict:
          homology_lengths.append(len(homology_dict[key].seq))
homology_lengths_len = len(homology_lengths)
print(homology_lengths_len)

if args.verbose:
    print('\033[0;34m' + "Loaded " + length_folder + '\033[1;37m')
    
# Open length restricted fasta
# Extract the number of reads contained within this file
length_dict = SeqIO.to_dict(SeqIO.parse(length_folder, "fasta"))

length_lengths = []
for key in length_dict:
    length_lengths.append(len(length_dict[key].seq))
length_lengths_len = len(length_dict)
print(length_lengths_len)

number_of_reads = pd.DataFrame([[full_lengths_len, homology_lengths_len, length_lengths_len]], columns=['Original # reads','# reads after homology filtering','# reads after length filtering'])
number_of_reads.to_csv(args.full_file[:9]+'Stats'+args.full_file[21:-12]+'read_numbers.csv', index=False)

if args.verbose:
    print('\033[0;32m' + "Number of reads file saved to " + args.full_file[:9]+'Stats'+args.full_file[21:-12]+'read_numbers.csv' + '\033[1;37m')

fields = [full_lengths_len, homology_lengths_len, length_lengths_len]
if args.full_file[40:-13] != 'unclassified':
    if '20171212_FAH18688/barcode10' not in args.full_file[22:-13] and '20171207_FAH18654/barcode10' not in args.full_file[22:-13]:
        with open(args.full_file[:9]+'Stats/read_numbers.csv','a+') as fd:
            writer = csv.writer(fd)
            writer.writerow(fields)
print('\033[0;35m'+'END'+'\033[1;37m')
