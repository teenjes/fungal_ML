
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
Enter the analysis directory to generate a dataframe that extracts the details
for each sample
""")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
parser.add_argument("input_directory", help="Run from honours folder and input analysis/")
args = parser.parse_args()

suspected = {'barcode01': 's__Puccinia_striiformis',
             'barcode02': 's__Zymoseptoria_tritici',
             'barcode03': 's__unknown', # healthy wheat
             'barcode04': 's__Puccinia_striiformis,s__Zymoseptoria_tritici', # both Pst + Zymo
             'barcode05': 's__Pyrenophora_tritici-repentis',
             'barcode06': 's__unknown', # healthy, resistant wheat, 2nd year,
             'barcode07': 's__unknown', # healthy, susceptible wheat, 2nd year
             'barcode08': 's__Puccinia_striiformis',
             'barcode09': 's__Puccinia_striiformis,s__Pyrenophora_tritici-repentis', # both Pst + Pyr
             'barcode10': 's__Puccinia_striiformis,s__Zymoseptoria_tritici', # both Pst + Zymo
             'barcode11': 's__Zymoseptoria_tritici,s__Puccinia_striiformis', # Zymo + Pst?
             'barcode12': 'none'# is just water
            }


summary = pd.DataFrame(data=None, columns = ['species 1', 'species 2','# raw reads','# reads after homology filtering','# reads after length filtering','path to raw reads','path to homology filtering','path to length filtering'])

import glob
path = "/media/MassStorage/tmp/TE/honours/analysis/Length_Filtered/wheat/*/"
path_names = glob.glob(path)
for path in path_names:
    if path[-13:-1] != 'unclassified' and path[-10:-1] != 'barcode12':
        key = path[-10:-1]
        if args.verbose:
            print('\033[0;34m' + "Opened barcode " + '\033[0;35m' + key + '\033[1;37m')
        if len(suspected[key].split(',')) == 1:
            species1 = suspected[key].split(',')[0].split('__')[1]
            species2 = None
        elif len(suspected[key].split(',')) == 2:
            species1 = suspected[key].split(',')[0].split('__')[1]
            species2 = suspected[key].split(',')[1].split('__')[1]
        
        if args.verbose:
            print('\033[1;36m' + "BEGIN RAW"'\033[1;37m')
        raw_path = "analysis/Concatenated/wheat/"+key+"/merged.fastq"
        raw_reads = SeqIO.to_dict(SeqIO.parse(raw_path, "fastq"))
        raw_count = 0
        for entry in raw_reads:
            raw_count += 1
            
        if args.verbose:
            print('\033[1;36m' + "BEGIN HOMOLOGY"+ '\033[1;37m')
        homology_path = "analysis/Python_Processing/wheat/"+key+"/combined_test.paf"
        homology_paf = pd.read_csv(homology_path, sep='\t', header=None, engine='python')
        homology_reads = {}
        homology_count = 0
        for entry in homology_paf[0].unique():
            homology_count += 1
            
        if args.verbose:
            print('\033[1;36m' + "BEGIN LENGTH" + '\033[1;37m')
        length_path = "analysis/Length_Filtered/wheat/"+key+"/length_restricted_reads.fasta"
        length_reads = SeqIO.to_dict(SeqIO.parse(length_path, "fasta"))
        length_count = 0
        for entry in length_reads:
            length_count += 1
            
        if args.verbose:
            print('\033[1;36m' + "BEGIN USE" + '\033[1;37m')
        use_path = length_path
        use_count = length_count


        add = pd.DataFrame([[species1, species2, raw_count, homology_count, length_count, use_count, raw_path, homology_path, length_path, use_path]], columns = ['species1', 'species2','# raw reads','# reads after homology filtering','# reads after length filtering','# for use', 'path to raw reads','path to homology filtering','path to length filtering', 'path for use'], index=[key])
        summary = summary.append([add])
summary = summary.sort_index(axis=0)
summary = summary[['species1', 'species2','# raw reads','# reads after homology filtering','# reads after length filtering','# for use', 'path to raw reads','path to homology filtering','path to length filtering', 'path for use']]
summary.to_csv("analysis/Stats/wheat_reference_dataframe.csv")
if args.verbose:
            print('\033[0;34m' + "Reference Dataframe saved to " + '\033[0;35m' + "analysis/Stats/wheat_reference_dataframe.csv" + '\033[1;37m')
