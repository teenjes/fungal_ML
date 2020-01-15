
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
import csv
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
args = parser.parse_args()

if args.verbose:
    print('\033[0;31m' + "Input file is " + args.input_file + '\033[1;37m')

keys = []
test_dict = SeqIO.to_dict(SeqIO.parse(args.input_file,"fasta"))
for key in test_dict:
    keys.append(key)
emboss_test = str(test_dict[keys[0]].seq)
emboss_new = emboss_test.replace('N','').replace('n','').replace('\n','')

scores_matrix = []
test50 = SeqIO.to_dict(SeqIO.parse("analysis/Alignment/20171103_FAH15473/barcode02/100_reads.fasta","fasta"))
for key in test50:
    alignments = pairwise2.align.globalxx(emboss_new, test50[key].seq, score_only=True)
    scores_matrix.append(int(alignments))

print(scores_matrix)
mean = np.mean(scores_matrix)
median=np.median(scores_matrix)
f=open((args.input_file[:47]+"scores.txt"),"a+")
f.write(">%s\t%s\t%s\t\n" % (keys[0], mean, median))

tmp = pd.DataFrame(scores_matrix)
tmp.to_csv(args.input_file[:47]+"%s.csv" % keys[0],index=False,header=False)
