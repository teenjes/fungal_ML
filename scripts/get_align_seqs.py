
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord as SR
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import AlignIO
from Bio.Align import AlignInfo
import numpy as np
import pandas as pd
from pandas import DataFrame as df
# import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import subprocess
import os
from shutil import copy
import random
import warnings
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
# from ipysankeywidget import SankeyWidget
warnings.filterwarnings("ignore")
import argparse


WHITE='\033[1;37m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'
PURPLE='\033[1;35m'

parser = argparse.ArgumentParser(description="""
This program extracts a specified number of reads from a fasta
file and saves them to a new fasta file
""")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
parser.add_argument("input_file", help="The input file for extraction")
parser.add_argument("num_reads", help="The number of reads to extract")
args = parser.parse_args()

if args.verbose:
    print('\033[0;31m' + "Input file is " + args.input_file + '\033[1;37m')
    print('\033[0;31m' + "The number of extracted reads is " + '\033[0;32m' + args.num_reads + '\033[1;37m')
    
tmp_dict = SeqIO.to_dict(SeqIO.parse(args.input_file,"fasta"))
new_dict = tmp_dict.copy()

ids=[]
for key in tmp_dict:
    ids.append(key)
    
keys_list = random.sample(ids,k=int(args.num_reads))

for key in new_dict:
    if key not in keys_list:
        del tmp_dict[key]
SeqIO.write(tmp_dict.values(),(args.input_file[:9]+'Alignment'+args.input_file[24:-29]+ args.num_reads + '_reads.fasta'),'fasta')

if args.verbose:
    print('\033[0;34m' + "Ids file saved to " + '\033[0;35m' + (args.input_file[:9]+'Alignment'+args.input_file[24:-29]+'for_consensus.fasta') + '\033[1;37m')
