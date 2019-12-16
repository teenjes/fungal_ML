
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
from mothur_py import Mothur
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
The goal of this program is to extract the read ids of reads that
contain both the forward and reverse primer as exact matches
""")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
parser.add_argument("input_file", help="The input file for extraction")
args = parser.parse_args()

if args.verbose:
    print('\033[0;31m' + "Input file is " + args.input_file + '\033[1;37m')

m = Mothur()
m.pcr.seqs(fasta=args.input_file,oligos=args.input_file[:9]+"ITS_primers.oligos",pdiffs=0,rdiffs=0)
pcr_dict = SeqIO.to_dict(SeqIO.parse(args.input_file[:-5]+"pcr.fasta","fasta"))
ids = []
for key in pcr_dict:
    ids.append(key)
with open(args.input_file[:-29]+'ids.txt','w') as handle:
    handle.writelines("%s\n" % name for name in ids)
    
if args.verbose:
    print('\033[0;34m' + "Ids file saved to " + '\033[0;35m' + (args.input_file[:-29]+'ids.txt') + '\033[1;37m')
    print('\033[0;32m' + ("The number of reads is %s" % len(ids)) + '\033[1;37m')
    
    
tmp_dict = SeqIO.to_dict(SeqIO.parse(args.input_file,"fasta"))
new_dict = tmp_dict.copy()
if len(ids) > 199:
    keys_list = random.sample(ids,k=200)
elif 100 < len(ids) < 200:
    keys_list = random.sample(ids,k=100)
else:
    print('\033[1;37m' + "LOW READS")
for key in new_dict:
    if key not in keys_list:
        del tmp_dict[key]
SeqIO.write(tmp_dict.values(),(args.input_file[:9]+'Consensus'+args.input_file[24:-29]+'for_consensus.fasta'),'fasta')

if args.verbose:
    print('\033[0;34m' + "Ids file saved to " + '\033[0;35m' + (args.input_file[:9]+'Consensus'+args.input_file[24:-29]+'for_consensus.fasta') + '\033[1;37m')
