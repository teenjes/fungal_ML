
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
This program alters the title of entries in the consensus fasta
""")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
parser.add_argument("input_file", help="The input file for extraction")
args = parser.parse_args()

true_names = {'20171103_FAH15473/barcode01': 'Puccinia_striiformis_tritici',
             '20171103_FAH15473/barcode02': 'Zymoseptoria_tritici',
             '20171103_FAH15473/barcode03': 'Pyrenophora_tritici-repentis',
             '20171103_FAH15473/barcode04': 'Fusarium_oxysporum',
             '20171103_FAH15473/barcode05': 'Tuber_brumale',
             '20171103_FAH15473/barcode06': 'Cortinarius_globuliformis',
             '20171103_FAH15473/barcode07': 'Aspergillus_niger',
             '20171103_FAH15473/barcode08': 'Clavispora_lusitaniae',
             '20171103_FAH15473/barcode09': 'Cryptococcus_neoformans',
             '20171103_FAH15473/barcode10': 'Penicillium_chrysogenum',
             '20171103_FAH15473/barcode11': 'Rhodotorula_mucilaginosa',
             '20171103_FAH15473/barcode12': 'Scedosporium_boydii',
             '20171207_FAH18654/barcode01': 'Blastobotrys_proliferans',
             '20171207_FAH18654/barcode02': 'Candida_zeylanoides',
             '20171207_FAH18654/barcode03': 'Galactomyces_geotrichum',
             '20171207_FAH18654/barcode04': 'Kodamaea_ohmeri',
             '20171207_FAH18654/barcode05': 'Meyerozyma_guillermondii',
             '20171207_FAH18654/barcode06': 'Wickerhamomyces_anomalus',
             '20171207_FAH18654/barcode07': 'Yamadazyma_mexicana',
             '20171207_FAH18654/barcode08': 'Yamadazyma_scolyti',
             '20171207_FAH18654/barcode09': 'Yarrowia_lipolytica',
             '20171207_FAH18654/barcode11': 'Zygoascus_hellenicus',
             '20171207_FAH18654/barcode12': 'Aspergillus_flavus',
             '20171212_FAH18688/barcode01': 'Cryptococcus_zero',
             '20171212_FAH18688/barcode02': 'Aspergillus_sp.',
             '20171212_FAH18688/barcode03': 'CCL067',
             '20171212_FAH18688/barcode04': 'Diaporthe_sp.',
             '20171212_FAH18688/barcode05': 'Tapesia_yallundae_CCL031',
             '20171212_FAH18688/barcode06': 'Tapesia_yallundae_CCL029',
             '20171212_FAH18688/barcode07': 'Dothiorella_vidmadera',
             '20171212_FAH18688/barcode08': 'Quambalaria_cyanescens',
             '20171212_FAH18688/barcode09': 'Entoleuca_sp.',
             '20171212_FAH18688/barcode11': 'CCL060',
             '20171212_FAH18688/barcode12': 'CCL068',
             '20180108_FAH18647/barcode01': 'Saccharomyces_cerevisiae',
             '20180108_FAH18647/barcode02': 'Cladophialophora_sp.',
             '20180108_FAH18647/barcode03': 'Candida_albicans',
             '20180108_FAH18647/barcode04': 'Candida_metapsilosis',
             '20180108_FAH18647/barcode05': 'Candida_orthopsilosis',
             '20180108_FAH18647/barcode06': 'Candida_parapsilosis',
             '20180108_FAH18647/barcode07': 'Cryptococcus_gattii',
             '20180108_FAH18647/barcode08': 'Geotrichum_candidum',
             '20180108_FAH18647/barcode09': 'Kluyveromyces_lactis',
             '20180108_FAH18647/barcode10': 'Kluyveromyces_marxianus',
             '20180108_FAH18647/barcode11': 'Pichia_kudriavzevii',
             '20180108_FAH18647/barcode12': 'Pichia_membranifaciens'}

if args.verbose:
    print("Input file is " + args.input_file + "\n")

with open(args.input_file) as original, open(args.input_file[:-25]+'final_custom_database_labelled.fasta', 'w+') as labelled:
    for line in original:
        if line.startswith(">"):
            barcode = line.split('/')[2]+'/'+line.split('/')[3]
            species = true_names[barcode]
            newline = '>'+species+'\n'
        else:
            newline = line
        labelled.write(newline)
original.close()
labelled.close()
