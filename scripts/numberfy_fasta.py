
import Bio
from Bio import SeqIO
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
This creates a machine-learning friendly version of each length restricted fasta file
""")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
parser.add_argument("input_file", help="The input file for extraction")
args = parser.parse_args()

if args.verbose:
    print("Input file is " + args.input_file + "\n")

fasta = SeqIO.to_dict(SeqIO.parse(args.input_file, "fasta"))

fasta_numbers = fasta.copy()

for key in fasta:
    seq = str(fasta[key].seq[30:-30]).replace("A",'0').replace("C",'1').replace("G",'2').replace("T",'3')
#     if len(seq) < max(total_lens):
#         seq = seq + '4'*(max(total_lens)-len(seq))
    fasta_numbers[key].seq = Seq(seq)

with open(args.input_file[:-6]+'_numbers.fasta', "w") as output_handle:
    SeqIO.write(fasta_numbers.values(), output_handle, 'fasta')
