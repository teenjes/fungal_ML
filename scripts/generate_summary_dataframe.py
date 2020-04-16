
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

phylogeny = {'20171103_FAH15473/barcode01': 'k__Fungi;p__Basidiomycota;c__Pucciniomycetes;o__Pucciniales;f__Pucciniaceae;g__Puccinia;s__striiformis-tritici',
             '20171103_FAH15473/barcode02': 'k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Capnodiales;f__Mycosphaerellaceae;g__Zymoseptoria;s__tritici',
             '20171103_FAH15473/barcode03': 'k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Pleosporales;f__Pleosporaceae;g__Pyrenophora;s__tritici-repentis',
             '20171103_FAH15473/barcode04': 'k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__oxysporum',
             '20171103_FAH15473/barcode05': 'k__Fungi;p__Ascomycota;c__Pezizomycetes;o__Pezizales;f__Tuberaceae;g__Tuber;s__brumale',
             '20171103_FAH15473/barcode06': 'k__Fungi;p__Basidiomycota;c__Agaricomyectes;o__Agaricales;f__Cortinariaceae;g__Cortinarius;s__globuliformis',
             '20171103_FAH15473/barcode07': 'k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Aspergillus;s__niger',
             '20171103_FAH15473/barcode08': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Metschnikowiaceae;g__Clavispora;s__lusitaniae',
        '20171103_FAH15473/barcode09': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Kluyveromyces;s__unidentified', # formerly Cryptococcus neoformans, matched in both UNITE and tree to Kluyveromyces
             '20171103_FAH15473/barcode10': 'k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Penicillium;s__chrysogenum',
             '20171103_FAH15473/barcode11': 'k__Fungi;p__Basidiomycota;c__Microbotryomycetes;o__Sporidiobolales;f__Sporidiobolaceae;g__Rhodotorula;s__mucilaginosa',
             '20171103_FAH15473/barcode12': 'k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Microascales;f__Microascaceae;g__Scedosporium;s__boydii',
             '20171207_FAH18654/barcode01': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Trichomonascaceae;g__Blastobotrys;s__proliferans', # matches Sugiyamella novakii in UNITE
             '20171207_FAH18654/barcode02': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Debaryomycetaceae;g__Debaryomyces;s__unidentified', # formerly Candida zeylanoides, matched in both UNITE and tree as Debaryomces/Meyerozyma
        '20171207_FAH18654/barcode03': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Dipodascaceae;g__Galactomyces;s__geotrichum', # tree puts it further away from other samples of its type, UNITE calls it Dipodascus
                                     # May have larger issues with Geotrichum candidum - maybe leave out?
             '20171207_FAH18654/barcode04': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Metschnikowiaceae;g__Kodamaea;s__ohmeri',
             '20171207_FAH18654/barcode05': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Debaryomycetaceae;g__Meyerozyma;s__guilliermondii',
             '20171207_FAH18654/barcode06': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Phaffomycetaceae;g__Wickerhamomyces;s__anomalus',
             '20171207_FAH18654/barcode07': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Debaryomycetaceae;g__Yamadazyma;s__mexicana',
             '20171207_FAH18654/barcode08': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Debaryomycetaceae;g__Yamadazyma;s__scolyti',
             '20171207_FAH18654/barcode09': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Dipodascaceae;g__Yarrowia;s__lipolytica',
             '20171207_FAH18654/barcode11': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Trichomonascaceae;g__Zygoascus;s__hellenicus', # matches Sugiyamella novakii in UNITE
             '20171207_FAH18654/barcode12': 'k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Aspergillus;s__flavus',
             
             '20171212_FAH18688/barcode01': 'k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Tremellaceae;g__Cryptococcus;s__zero',
             '20171212_FAH18688/barcode02': 'k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Aspergillus;s__unidentified',
         '20171212_FAH18688/barcode03': 'k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Diaporthales;f__Diaporthaceae;g__Diaporthe;s__CCL067',
             '20171212_FAH18688/barcode04': 'k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Diaporthales;f__Diaporthaceae;g__Diaporthe;s__unidentified',
             '20171212_FAH18688/barcode05': 'k__Fungi;p__Ascomycota;c__Leotiomycetes;o__Helotiales;f__Dermateaceae;g__Oculimacula;s__yallundae-CCL031',
             '20171212_FAH18688/barcode06': 'k__Fungi;p__Ascomycota;c__Leotiomycetes;o__Helotiales;f__Dermateaceae;g__Oculimacula;s__yallundae-CCL029',
             '20171212_FAH18688/barcode07': 'k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Botryosphaeriales;f__Botryosphaeriaceae;g__Dothiorella;s__vidmadera',
             '20171212_FAH18688/barcode08': 'k__Fungi;p__Basidiomycota;c__Exobasidiomycetes;o__Microstromatales;f__Quambalariaceae;g__Quambalaria;s__cyanescens',
             '20171212_FAH18688/barcode09': 'k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Xylariales;f__Xylariaceae;g__Entoleuca;s__unidentified',
                 '20171212_FAH18688/barcode11': 'k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Diaporthales;f__Gnomoniaceae;g__Asteroma;s__CCL060', # may be asteroma or Diaporthe
                 '20171212_FAH18688/barcode12': 'k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Diaporthales;f__Gnomoniaceae;g__Asteroma;s__CCL068',
             '20180108_FAH18647/barcode01': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Saccharomyces;s__cerevisiae',
             '20180108_FAH18647/barcode02': 'k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Chaetothyriales;f__Herpotrichiellaceae;g__Cladophialophora;s__unidentified',
             '20180108_FAH18647/barcode03': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Candida;s__albicans',
             '20180108_FAH18647/barcode04': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Candida;s__metapsilosis',
             '20180108_FAH18647/barcode05': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Candida;s__orthopsilosis',
             '20180108_FAH18647/barcode06': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Candida;s__parapsilosis',
#          '20180108_FAH18647/barcode07': 'k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Tremellaceae;g__Cryptococcus;s__gattii', # not including
#              '20180108_FAH18647/barcode08': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Dipodascaceae;g__Geotrichum;s__candidum', # not including
             '20180108_FAH18647/barcode09': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Candida;s__unidentified', # formerly Kluyveromyces lactis
             '20180108_FAH18647/barcode10': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Kluyveromyces;s__marxianus',
             '20180108_FAH18647/barcode11': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Pichiaceae;g__Pichia;s__kudriavzevii',
             '20180108_FAH18647/barcode12': 'k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Pichiaceae;g__Pichia;s__membranifaciens'}

summary = pd.DataFrame(data=None, columns = ['species','genus','family','order','class','phylum','kingdom','# raw reads','# reads after homology filtering','# reads after length filtering','path to raw reads','path to homology filtering','path to length filtering'])
import glob
path = "/media/MassStorage/tmp/TE/honours/analysis/Consensus/*/*/"
path_names = glob.glob(path)
for path in path_names:
    if path[-13:-1] != 'unclassified' and path[-28:-1] != '20171207_FAH18654/barcode10' and path[-28:-1] != '20171212_FAH18688/barcode10' and path[-28:-1] != '20180108_FAH18647/barcode07' and path[-28:-1] != '20180108_FAH18647/barcode08':
        key = path[-28:-1]
        if args.verbose:
            print('\033[0;34m' + "Opened barcode " + '\033[0;35m' + key + '\033[1;37m')
        species_ = phylogeny[key].split(';')[6].split('__')[1].lower()
        genus_ = phylogeny[key].split(';')[5].split('__')[1].lower()
        family_ = phylogeny[key].split(';')[4].split('__')[1].lower()
        order_ = phylogeny[key].split(';')[3].split('__')[1].lower()
        class_ = phylogeny[key].split(';')[2].split('__')[1].lower()
        phylum_ = phylogeny[key].split(';')[1].split('__')[1].lower()
        kingdom_ = phylogeny[key].split(';')[0].split('__')[1].lower()
        
        if args.verbose:
            print('\033[1;36m' + "BEGIN RAW"'\033[1;37m')
        raw_path = "analysis/Concatenated/"+key+"/merged.fastq"
        raw_reads = SeqIO.to_dict(SeqIO.parse(raw_path, "fastq"))
        raw_count = 0
        for entry in raw_reads:
            raw_count += 1
            
        if args.verbose:
            print('\033[1;36m' + "BEGIN HOMOLOGY"+ '\033[1;37m')
        homology_path = "analysis/Python_Processing/"+key+"/combined_test.paf"
        homology_paf = pd.read_csv(homology_path, sep='\t', header=None, engine='python')
        homology_reads = {}
        homology_count = 0
        for entry in homology_paf[0].unique():
            homology_count += 1
            
        if args.verbose:
            print('\033[1;36m' + "BEGIN LENGTH" + '\033[1;37m')
        length_path = "analysis/Length_Filtered/"+key+"/length_restricted_reads.fasta"
        length_reads = SeqIO.to_dict(SeqIO.parse(length_path, "fasta"))
        length_count = 0
        for entry in length_reads:
            length_count += 1
            
        if args.verbose:
            print('\033[1;36m' + "BEGIN USE" + '\033[1;37m')
        use_path = "analysis/Length_Filtered/"+key+"/length_restricted_for_use.fasta"
        use_reads = SeqIO.to_dict(SeqIO.parse(use_path, "fasta"))
        use_count = 0
        for entry in use_reads:
            use_count += 1


        add = pd.DataFrame([[species_, genus_, family_, order_, class_, phylum_, kingdom_, raw_count, homology_count, length_count, use_count, raw_path, homology_path, length_path, use_path]], columns = ['species','genus','family','order','class','phylum','kingdom','# raw reads','# reads after homology filtering','# reads after length filtering','# for use', 'path to raw reads','path to homology filtering','path to length filtering', 'path for use'], index=[key])
        summary = summary.append([add])
summary = summary.sort_index(axis=0)
summary = summary[['species','genus','family','order','class','phylum','kingdom','# raw reads','# reads after homology filtering','# reads after length filtering','# for use', 'path to raw reads','path to homology filtering','path to length filtering', 'path for use']]
summary.to_csv("analysis/Stats/reference_dataframe.csv")
if args.verbose:
            print('\033[0;34m' + "Reference Dataframe saved to " + '\033[0;35m' + "analysis/Stats/reference_dataframe.csv" + '\033[1;37m')
