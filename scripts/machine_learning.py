
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

def max_seq_len(SeqIO_dict):
    """
    Function takes a SeqIO_dict and returns the lengths of the
    longest sequence
    """
    total_lens = []
    for key in SeqIO_dict.keys():
        total_lens.append(len(SeqIO_dict[key].seq))
    return max(total_lens)

def numberfy(SeqIO_dict, seq_len, nsubsample):
    """
    Take SeqIO_dict and return SeqIO_dict were bases have been replaced
    with numbers
    ACGT- replaced with 01234
    Take the seq_len each sequence should have
    """
    num_dict = {}
    
    keys = list(SeqIO_dict.keys())
    randkeys = random.sample(keys, k=nsubsample)
    
    
    for key in randkeys:
        seq = str(SeqIO_dict[key].seq).replace("A",'0 ')\
        .replace("C",'1 ').replace("G",'2 ').replace("T",'3 ')
        seq_new = seq + '4 '*(seq_len -int(len(seq)/2))
        num_dict[key] = list(map(int, seq_new.split(' ')[:-1]))
    return num_dict

parser = argparse.ArgumentParser(description="""
Script aims to read in the reference_dataframe file, select a taxonomic
level and group, and read the path to the location of that data. It then
prepares data for machine learning by converting base pair coding to numerical
encoding, pads it out and then runs the algorithm
""")
parser.add_argument("ref_df_fn", help="File path to the reference dataframe")
parser.add_argument("data_root", help="Root folder for analysis/")
parser.add_argument("--tax_rank", "-r", help="taxonomic rank for analysis")
parser.add_argument("--name", "-n", help="name of rank to select from")
parser.add_argument("--n_reads", "-c", help="count of reads per class")
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
    
# assign flagged variables as lower case
tax_rank = args.tax_rank.lower()
name = args.name.lower()

if args.verbose:
    print('\033[1;34m' + "Reference dataframe is at " + ref_df_fn + '\033[0m')
    print('\033[1;34m' + "Root directory is at " + data_root + '\033[0m')
    print('\033[1;34m' + "Tax Rank is " + tax_rank + '\033[0m')
    print('\033[1;34m' + "Name is " + name + '\033[0m')
    print('\033[1;34m' + "Count of reads per sample is", n_reads,'\033[0m')

# read in the reference dataframe from the argument path
ref_df = pd.read_csv(ref_df_fn, index_col=None)

# check whether the reference dataframe implies there are enough reads
# to continue given n_reads
try:
    if ref_df[ref_df["# reads after length filtering"] \
              < n_reads].shape[0] > 0 :
        print("These species need more reads.")
        print(ref_df[ref_df["# reads after length filtering"] \
              < n_reads])
        #exit()
except:
    print('Check %s to have the wanted column names' % ref_df_fn)

# assign the indices of the reference_dataframe as keys to a dictionary
# where the values are that index's path's dataframe
indices = ref_df[ref_df[tax_rank] == name].index
SeqIO_dicts = {}
for index in indices:
    fasta_path = ref_df.loc[index, 'path to length filtering']
    try:
        SeqIO_dicts[index] = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    except:
        print('Check location of fasta files')
        print(fasta_path, "does not exist")
        
# determine the maximum sequence length of accepted sequences
total_lens = []
for key, value in SeqIO_dicts.items():
    total_lens.append(max_seq_len(value))
print('\033[0;32m'+"The maximum sequence length of all sampled sequences is"+ '\033[1;37m',max(total_lens),'\033[0m')

# randomly subsample n_reads number of reads from each index's corresponding
# set of reds, convert base pair coding to numerical coding and 
# pad to the max sequence length
numSeqIO_dicts = {}
max_len = max(total_lens)
for key, value in SeqIO_dicts.items():
    numSeqIO_dicts[key] = numberfy(value, max_len, n_reads)

# append the numberfy'd sequences to a numpy array
seq_list = []
for index in indices:
    seq_list.append(np.array(list(numSeqIO_dicts[index].values())))
seq_comb = np.concatenate(seq_list, axis = 0)

# determine the number of classes and generate an array of ids
num_class = len(numSeqIO_dicts.keys())
ids_comb = np.zeros( (n_reads*num_class,num_class) )
for i in range(0, num_class):
    ids_comb[i*n_reads:(i+1)*n_reads,i] = 1

print(ids_comb)
