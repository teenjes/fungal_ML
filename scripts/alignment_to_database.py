
"""
Script aims to read in the reference_dataframe file, select a taxonomic
level and group, and read the path to the location of that data. It then
prepares data for machine learning by converting base pair coding to numerical
encoding, pads it out and then runs the algorithm
"""

import pandas as pd
from Bio import SeqIO
import numpy as np
import os
import random
import argparse
import math
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(description="""
Script aims to read in the reference_dataframe file, and for each
sample, randomly subsample n_reads number of reads and save the
file, and save the keys used as as separate file. The script ends there.

Outside the script these files will be utilised in minimap2 to
map against a specified database, produce an output paf file
and then remove the file to conserve space. 
""")
parser.add_argument("ref_df_fn", help="File path to the reference dataframe")
parser.add_argument("data_root", help="Root folder for analysis/")
parser.add_argument("--n_reads", "-c", help="count of reads per class")
parser.add_argument("--d_type", "-d", help="type of database - CUSTOM or UNITE")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
args = parser.parse_args()

# assign required arguments to variables
ref_df_fn = args.ref_df_fn
data_root = args.data_root

# assign a number of reads per class
n_reads = int(args.n_reads)

# test to make sure both required file paths are input
try:
    os.path.exists(ref_df_fn)
except:
    print('Cannot find %s' % ref_df_fn)
try:
    os.path.exists(data_root)
except:
    print('Cannot find %s' % data_root)
    
if args.verbose:
    print('\033[1;34m' + "Reference dataframe is at " + ref_df_fn + '\033[0m')
    print('\033[1;34m' + "Root directory is at " + data_root + '\033[0m')
    print('\033[1;34m' + "Count of reads per sample is", n_reads,'\033[0m')
    
# read in the reference dataframe from the argument path
ref_df = pd.read_csv(ref_df_fn, index_col=None)

# check whether the reference dataframe implies there are enough reads
# to continue given n_reads
try:
    if ref_df[ref_df["# for use"] \
              < n_reads].shape[0] > 0 :
        print("These species need more reads.")
        print(ref_df[ref_df["# for use"] \
              < n_reads])
        exit()
except:
    print('Check %s to have the wanted column names' % ref_df_fn)
   

for index, row in ref_df.iterrows():
    fasta_dict = SeqIO.to_dict(SeqIO.parse(row['path for use'], "fasta"))
    key_list = []
    for key in fasta_dict:
        key_list.append(key)
    print(len(key_list))
    keys_list = random.sample(key_list,k=n_reads)
    new_dict = {}
    for key in keys_list:
        new_dict[key] = fasta_dict[key]
    if args.d_type.upper() == 'CUSTOM':
        SeqIO.write(new_dict.values(), data_root+"database_mapping/custom/%s_%s.fasta" % (row['genus'],row['species']), "fasta")
        with open(data_root+'database_mapping/custom/%s_%s_keys.csv'% (row['genus'],row['species']), 'w+') as f:
            for key in keys_list:
                f.write("%s\n"%(key))
    elif args.d_type.upper() == 'UNITE':
        SeqIO.write(new_dict.values(), "analysis/database_mapping/unite/%s_%s.fasta" % (row['genus'],row['species']), "fasta")
        with open(data_root+'database_mapping/unite/%s_%s_keys.csv'% (row['genus'],row['species']), 'w+') as f:
            for key in keys_list:
                f.write("%s\n"%(key))
