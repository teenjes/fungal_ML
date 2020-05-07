
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

suspected = {'barcode01': 'k__Fungi;p__Basidiomycota;c__Pucciniomycetes;o__Pucciniales;f__Pucciniaceae;g__Puccinia;s__Puccinia_striiformis',
             'barcode02': 'k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Capnodiales;f__Mycosphaerellaceae;g__Zymoseptoria;s__Zymoseptoria_tritici',
             'barcode05': 'k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Pleosporales;f__Pleosporaceae;g__Pyrenophora;s__Pyrenophora_tritici-repentis',
            }


summary = pd.DataFrame(data=None, columns = ['species','genus','family','order','class','phylum','kingdom','# raw reads','# reads after homology filtering','# reads after length filtering','path to raw reads','path to homology filtering','path to length filtering'])

import glob
path = "/media/MassStorage/tmp/TE/honours/analysis/Length_Filtered/wheat/*/"
path_names = glob.glob(path)
for path in path_names:
    if path[-10:-1] == 'barcode01' or path[-10:-1] == 'barcode02' or path[-10:-1] == 'barcode05':
        key = path[-10:-1]
        if args.verbose:
            print('\033[0;34m' + "Opened barcode " + '\033[0;35m' + key + '\033[1;37m')
        species_ = suspected[key].split(';')[6].split('__')[1].lower()
        genus_ = suspected[key].split(';')[5].split('__')[1].lower()
        family_ = suspected[key].split(';')[4].split('__')[1].lower()
        order_ = suspected[key].split(';')[3].split('__')[1].lower()
        class_ = suspected[key].split(';')[2].split('__')[1].lower()
        phylum_ = suspected[key].split(';')[1].split('__')[1].lower()
        kingdom_ = suspected[key].split(';')[0].split('__')[1].lower()
        
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


        add = pd.DataFrame([[species_, genus_, family_, order_, class_, phylum_, kingdom_, raw_count, homology_count, length_count, use_count, raw_path, homology_path, length_path, use_path]], columns = ['species','genus','family','order','class','phylum','kingdom', '# raw reads','# reads after homology filtering','# reads after length filtering','# for use', 'path to raw reads','path to homology filtering','path to length filtering', 'path for use'], index=[key])
        summary = summary.append([add])
summary = summary.sort_index(axis=0)
summary = summary[['species','genus','family','order','class','phylum','kingdom','# raw reads','# reads after homology filtering','# reads after length filtering','# for use', 'path to raw reads','path to homology filtering','path to length filtering', 'path for use']]
summary.to_csv("analysis/Stats/wheat_reference_dataframe.csv")
if args.verbose:
            print('\033[0;34m' + "Reference Dataframe saved to " + '\033[0;35m' + "analysis/Stats/wheat_reference_dataframe.csv" + '\033[1;37m')
