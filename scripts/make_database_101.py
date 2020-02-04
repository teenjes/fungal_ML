
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
This program writes a cleaned consensus to a database file
""")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
parser.add_argument("input_file", help="The input file for extraction")
args = parser.parse_args()

if args.verbose:
    print("Input file is " + args.input_file + "\n")

keys = []
test_dict = SeqIO.to_dict(SeqIO.parse(args.input_file,"fasta"))
for key in test_dict:
    keys.append(key)
emboss_test = str(test_dict[keys[0]].seq)
emboss_new = emboss_test.replace('N','').replace('n','').replace('\n','')
    
with open("database/custom_database_101.fasta","a") as handle:
    handle.write(">%s\n" % keys[0] + emboss_new + '\n')
