
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
    print("Input file is " + args.input_file + "\n")

f=open(args.input_file,"r")
if f.mode == "r":
    contents=f.read()

tmp=contents.replace("\t",",").split('\n')
tmp_dict = {}
for line in tmp[:-1]:
    tmp_dict[line.split(",")[0]] = str(line.split(",")[5])

count_dict = {}
correct = 0
incorrect = 0
for item in tmp_dict:
    if tmp_dict[item] not in count_dict:
        count_dict[tmp_dict[item]] = 1
    else:
        count_dict[tmp_dict[item]] = count_dict[tmp_dict[item]] + 1
        
tmp = pd.DataFrame.from_dict(count_dict,orient='index',columns=["Count"])
tmp["Percentage Match"] = tmp.apply(lambda row: 100*row.Count/1000,axis=1)
tmp.index.names = ['analysis/Consensus/'+args.input_file[19:-18]]
tmp = tmp.sort_values(by="Count",ascending=False)

if args.verbose:
    print(tmp)

tmp.to_csv(args.input_file[:-18]+'match_distribution.csv',sep=',')
