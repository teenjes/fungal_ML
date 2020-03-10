
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
Input the sequencing_summary.txt file from Basecalled
Read in the length_restricted_reads.fasta from Length_Filtered
Save to Consensus
""")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
parser.add_argument("input_file", help="The input file for extraction")
args = parser.parse_args()

if args.verbose:
    print('\033[0;31m' + "Input file is " + args.input_file + '\033[1;37m')

    
tmp = pd.read_csv(args.input_file,sep='\t', header=None, names=['filename', 'read_id', 'run_id', 'batch_id', 'channel', 'mux', 'start_time', 'duration', 'num_events', 'passes_filtering', 'template_start', 'num_events_template', 'template_duration', 'sequence_length_template', 'mean_qscore_template', 'strand_score_template', 'median_template', 'mad_template'], engine='python')
tmpd = SeqIO.to_dict(SeqIO.parse("Length_Filtered"+args.input_file[10:-22]+"length_restricted_reads.fasta", "fasta"))
keys = []
for key in tmpd:
    keys.append(key)
id_score = pd.DataFrame(data=None,columns=['read_id','mean_qscore_template'])
tmp_keys = {}
for i in range(1,len(tmp)):
    tmp_keys[tmp['read_id'][i]] = tmp['mean_qscore_template'][i]
for key in keys:
    if key in tmp_keys:
        tmpf = pd.DataFrame(data=[[key,tmp_keys[key]]],columns=['read_id','mean_qscore_template'])
    id_score = id_score.append(tmpf)
id_score = id_score.sort_values(by=['mean_qscore_template'])
new_ids_list = id_score['read_id'][:100]

new_dict = {}
for key in new_ids_list:
    new_dict[key] = tmpd[key]
SeqIO.write(new_dict.values(), "Consensus"+args.input_file[10:-22]+"high_qscore_100.fasta", "fasta")
    

if args.verbose:
    print('\033[0;34m' + "Reads saved to " + '\033[0;35m' + "Consensus"+args.input_file[10:-22]+"high_qscore_100.fasta" + '\033[1;37m')
