
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
parser.add_argument("input_file", help="The input file for extraction")
args = parser.parse_args()

phylogeny = {'20171103_FAH15473/barcode01': 'Basidiomycota/Pucciniomycetes/Pucciniales/Pucciniaceae/Puccinia/striiformis_tritici',
             '20171103_FAH15473/barcode02': 'Ascomycota/Dothideomycetes/Capnodiales/Mycosphaerellaceae/Zymoseptoria/tritici',
             '20171103_FAH15473/barcode03': 'Ascomycota/Dothideomycetes/Pleosporales/Pleospraceae/Pyrenophora/tritici-repentis',
             '20171103_FAH15473/barcode04': 'Ascomycota/Sordariomycetes/Hypocreales/Nectriaceae/Fusarium/oxysporum',
             '20171103_FAH15473/barcode05': 'Ascomycota/Pezizomycetes/Pezizales/Tuberaceae/Tuber/brumale',
             '20171103_FAH15473/barcode06': 'Basidiomycota/Agaricomycetes/Agaricales/Cortinariaceae/Cortinarius/globuliformis',
             '20171103_FAH15473/barcode07': 'Ascomycota/Eurotiomycetes/Eurotiales/Trichocomaceae/Aspergillus/niger',
             '20171103_FAH15473/barcode08': 'Ascomycota/Saccharomycetes/Saccharomycetales/Metschnikowiaceae/Clavispora/lusitaniae',
                 '20171103_FAH15473/barcode09': 'Basidiomycota/Tremellomycetes/Tremellales/Tremellaceae/Cryptococcus/neoformans',
             '20171103_FAH15473/barcode10': 'Ascomycota/Eurotiomycetes/Eurotiales/Trichocomaceae/Penicillium/chrysogenum',
             '20171103_FAH15473/barcode11': 'Basidiomycota/Microbotryomycetes/Sporidiobolales/Sporidiobolaceae/Rhodotorula/mucilaginosa',
             '20171103_FAH15473/barcode12': 'Ascomycota/Sordariomycetes/Microascales/Microascaceae/Scedosporium/boydii',
             '20171207_FAH18654/barcode01': 'Ascomycota/Saccharomycetes/Saccharomycetales/Trichomonascaceae/Blastobotrys/proliferans',
                 '20171207_FAH18654/barcode02': 'Ascomycota/Saccharomycetes/Saccharomycetales/Saccharomycetaceae/Candida/zeylanoides',
                 '20171207_FAH18654/barcode03': 'Ascomycota/Saccharomycetes/Saccharomycetales/Dipodascaceae/Galactomyces/geotrichum',
             '20171207_FAH18654/barcode04': 'Ascomycota/Saccharomycetes/Saccharomycetales/Metschnikowiaceae/Kodamaea/ohmeri',
             '20171207_FAH18654/barcode05': 'Ascomycota/Saccharomycetes/Saccharomycetales/Debaryomycetaceae/Meyerozyma/guillermondii',
             '20171207_FAH18654/barcode06': 'Ascomycota/Saccharomycetes/Saccharomycetales/Phaffomycetaceae/Wickerhamomyces/anomalus',
             '20171207_FAH18654/barcode07': 'Ascomycota/Saccharomycetes/Saccharomycetales/Debaryomycetaceae/Yamadazyma/mexicana',
             '20171207_FAH18654/barcode08': 'Ascomycota/Saccharomycetes/Saccharomycetales/Debaryomycetaceae/Yamadazyma/scolyti',
             '20171207_FAH18654/barcode09': 'Ascomycota/Saccharomycetes/Saccharomycetales/Dipodascaceae/Yarrowia/lipolytica',
             '20171207_FAH18654/barcode11': 'Ascomycota/Saccharomycetes/Saccharomycetales/Trichomonascaceae/Zygoascus/hellenicus',
             '20171207_FAH18654/barcode12': 'Ascomycota/Eurotiomycetes/Eurotiales/Trichocomaceae/Aspergillus/flavus',
             '20171212_FAH18688/barcode01': 'Basidiomycota/Tremellomycetes/Tremellales/Tremellaceae/Cryptococcus/zero',
             '20171212_FAH18688/barcode02': 'Ascomycota/Eurotiomycetes/Eurotiales/Trichocomaceae/Aspergillus/sp.',
                 '20171212_FAH18688/barcode03': 'Ascomycota/Sordariomycetes/Diaporthales/Diaporthaceae/Diaporthe/sp._CCL067',
             '20171212_FAH18688/barcode04': 'Ascomycota/Sordariomycetes/Diaporthales/Diaporthaceae/Diaporthe/sp.',
             '20171212_FAH18688/barcode05': 'Ascomycota/Leotiomycetes/Helotiales/Dermateaceae/Oculimacula/yallundae_CCL031',
             '20171212_FAH18688/barcode06': 'Ascomycota/Leotiomycetes/Helotiales/Dermateaceae/Oculimacula/yallundae_CCL029',
             '20171212_FAH18688/barcode07': 'Ascomycota/Dothideomycetes/Botryosphaeriales/Botryosphaeriaceae/Dothiorella/vidmadera',
             '20171212_FAH18688/barcode08': 'Exobasidiomycetes/Microstromatales/Quambalariaceae/Quambalaria/cyanescens',
             '20171212_FAH18688/barcode09': 'Ascomycota/Sordariomycetes/Xylariales.Xylariaceae/Entoleuca/sp.',
                 '20171212_FAH18688/barcode11': 'Ascomycota/Sordariomycetes/Diaporthales/Diaporthaceae/Diaporthe/sp._CCL060',
                 '20171212_FAH18688/barcode12': 'Ascomycota/Sordariomycetes/Diaporthales/Diaporthaceae/Diaporthe/sp._CCL068',
             '20180108_FAH18647/barcode01': 'Ascomycota/Saccharomycetes/Saccharomycetales/Saccharomycetaceae/Saccharomyces/cerevisiae',
             '20180108_FAH18647/barcode02': 'Ascomycota/Eurotiomycetes/Chaetothyriales/Herpotrichiellaceae/Cladophialophora/sp.',
             '20180108_FAH18647/barcode03': 'Ascomycota/Saccharomycetes/Saccharomycetales/Saccharomycetaceae/Candida/albicans',
             '20180108_FAH18647/barcode04': 'Ascomycota/Saccharomycetes/Saccharomycetales/Saccharomycetaceae/Candida/metapsilosis',
             '20180108_FAH18647/barcode05': 'Ascomycota/Saccharomycetes/Saccharomycetales/Saccharomycetaceae/Candida/orthopsilosis',
             '20180108_FAH18647/barcode06': 'Ascomycota/Saccharomycetes/Saccharomycetales/Saccharomycetaceae/Candida/parapsilosis',
                 '20180108_FAH18647/barcode07': 'Basidiomycota/Tremellomycetes/Tremellales/Tremellaceae/Cryptococcus/gattii',
             '20180108_FAH18647/barcode08': 'Ascomycota/Saccharomycetes/Saccharomycetales/Dipodascaceae/Geotrichum/candidum',
             '20180108_FAH18647/barcode09': 'Ascomycota/Saccharomycetes/Saccharomycetales/Saccharomycetaceae/Kluyveromyces/lactis',
             '20180108_FAH18647/barcode10': 'Ascomycota/Saccharomycetes/Saccharomycetales/Saccharomycetaceae/Kluyveromyces/marxianus',
             '20180108_FAH18647/barcode11': 'Ascomycota/Saccharomycetes/Saccharomycetales/Pichiaceae/Pichia/kudriavzevii',
             '20180108_FAH18647/barcode12': 'Ascomycota/Saccharomycetes/Saccharomycetales/Pichiaceae/Pichia/membranifaciens'}
input_file = args.input_file

print(input_file+'\n')

f=open(input_file,"r")
if f.mode == "r":
    contents=f.read()  
tmp=contents.replace("\t",",").split('\n')
tmp
tmp_dict = {}
counter=0

for line in tmp[:-1]:
    tmp_dict[counter] = str(line.split(",")[5].split('_')[0])
    counter+=1
tmp2 = pd.DataFrame.from_dict(data=tmp_dict, orient='index')
counts = tmp2[0].value_counts()
tmp3 = counts.to_frame().reset_index(drop=False)

g=open("database/sh_taxonomy_qiime_ver8_dynamic_02.02.2019.txt", "r")
if g.mode == "r":
    contents2=g.read()
tmp4=contents2.replace("\t",",").split("\n")

tmp5 = []
for i in range(0,len(tmp4)):
    tmp5.append(tmp4[i].split(','))

tmp6 = pd.DataFrame(tmp5)

tmp_dict2 = {}
for i in range(0,len(tmp6[0])):
    tmp_dict2[tmp6[0][i].split("_")[0]]=tmp6[1][i]

    
tobewritten = pd.DataFrame(data=[],columns=['Count','Species Match','Genus Match','Family Match','Order Match','Class Match','Phylum Match'])

for a,b in tmp3[:7].itertuples(index=False):
    # a == id, b == count
    if a in tmp_dict2:
        tmp_specs = []
        tmp_specs.append(b)
        if len(tmp_dict2[a].split('s__')[-1].split('_')) > 1:
            tmp_specs.append(tmp_dict2[a].split('s__')[-1].split('_')[1])
            tmp_specs.append(tmp_dict2[a].split('s__')[-1].split('_')[0])
        else:
            tmp_specs.append(tmp_dict2[a].split('s__')[-1].split('_')[0])
            tmp_specs.append(tmp_dict2[a].split('s__')[-1].split('_')[0])
        tmp_specs.append(tmp_dict2[a].split('s__')[-2].split('g__')[-2].split('f__')[-1].split(';')[0])
        tmp_specs.append(tmp_dict2[a].split('s__')[-2].split('g__')[-2].split('f__')[-2].split('o__')[-1].split(';')[0])
        tmp_specs.append(tmp_dict2[a].split('s__')[-2].split('g__')[-2].split('f__')[-2].split('o__')[-2].split('c__')[-1].split(';')[0])
        tmp_specs.append(tmp_dict2[a].split('s__')[-2].split('g__')[-2].split('f__')[-2].split('o__')[-2].split('c__')[-2].split('p__')[-1].split(';')[0])
        tmpf = pd.DataFrame([tmp_specs],columns=['Count','Species Match','Genus Match','Family Match','Order Match','Class Match','Phylum Match'],index=[a])
        tobewritten = tobewritten.append(tmpf)

        
pd.set_option('display.max_rows', None)
print('\033[0;34m')
print(tobewritten)
print('\033[1;37m')       
# correct_list = ['N/A', phylogeny[input_file[19:-33]].split('/')[-1].split('_')[0], phylogeny[input_file[19:-33]].split('/')[-2],phylogeny[input_file[19:-33]].split('/')[-3],phylogeny[input_file[19:-33]].split('/')[-4],phylogeny[input_file[19:-33]].split('/')[-5],phylogeny[input_file[19:-33]].split('/')[-6]]
# correct = pd.DataFrame([correct_list], columns=['Count','Species Match','Genus Match','Family Match','Order Match','Class Match','Phylum Match'],index=['Suggested'])

# tps1 = correct[['Species Match']]
# tps2 = tobewritten[['Species Match']]
# tpg1 = correct[['Genus Match']]
# tpg2 = tobewritten[['Genus Match']]
# print(pd.merge(tps1,tps2,on='Species Match'))
# if pd.merge(tps1,tps2,on='Species Match').any().any():
#     print('\033[0;34m')
#     print(pd.merge(tps1,tps2,on='Species Match').any())
# else:
#     print('\033[0;31m')
#     print(pd.merge(tps1,tps2,on='Species Match').any())
# print('\033[1;37m')

# print(pd.merge(tpg1,tpg2,on='Genus Match'))
# if pd.merge(tpg1,tpg2,on='Genus Match').any().any():
#     print('\033[0;34m')
#     print(pd.merge(tpg1,tpg2,on='Genus Match').any())
# else:
#     print('\033[0;31m')
#     print(pd.merge(tpg1,tpg2,on='Genus Match').any())
# print('\033[1;37m')



# consensus = [0]*7
# count = 0
# total = 0
# match = 0
# for line in tobewritten:
#     if line != 'Count':
#         tmps = tobewritten[line].value_counts().to_frame()
#         if tmps.index[0] in correct_list:
#             consensus[count] = tmps.index[0]
#             print(consensus[count])
#         else:
#             print("NO")
#         count += 1
#     else:
#         consensus[count] = 'N/A'
#         count += 1

        
# correct = pd.DataFrame([correct_list], columns=['Count','Species Match','Genus Match','Family Match','Order Match','Class Match','Phylum Match'],index=['Suggested'])
# tobewritten = tobewritten.append(correct)

# consensus = pd.DataFrame([consensus], columns=['Count','Species Match','Genus Match','Family Match','Order Match','Class Match','Phylum Match'],index=['Consensus Match'])
# tobewritten = tobewritten.append([consensus])

# if (correct.values[0][1:2] == consensus.values[0][1:2]).all():
#     print('\033[0;31m'+"SPECIES MATCH"+'\033[1;37m')
# else:
#     print('\033[0;34m'+"NO SPECIES MATCH"+'\033[1;37m')
    
# if (correct.values[0][2:3] == consensus.values[0][2:3]).all():
#     print('\033[0;31m'+"GENUS MATCH"+'\033[1;37m')
# else:
#     print('\033[0;34m'+"NO GENUS MATCH"+'\033[1;37m')
    
# if (correct.values[0][3:4] == consensus.values[0][3:4]).all():
#     print('\033[0;31m'+"FAMILY MATCH"+'\033[1;37m')
# else:
#     print('\033[0;34m'+"NO FAMILY MATCH"+'\033[1;37m')
    
# if (correct.values[0][4:5] == consensus.values[0][4:5]).all():
#     print('\033[0;31m'+"ORDER MATCH"+'\033[1;37m')
# else:
#     print('\033[0;34m'+"NO ORDER MATCH"+'\033[1;37m')
    
# if (correct.values[0][5:6] == consensus.values[0][5:6]).all():
#     print('\033[0;31m'+"CLASS MATCH"+'\033[1;37m')
# else:
#     print('\033[0;34m'+"NO CLASS MATCH"+'\033[1;37m')
    
# if (correct.values[0][6:] == consensus.values[0][6:]).all():
#     print('\033[0;31m'+"PHYLUM MATCH"+'\033[1;37m')
# else:
#     print('\033[0;34m'+"NO PHYLUM MATCH"+'\033[1;37m')
    
# print('\n\n')


# tobewritten.to_csv(input_file[:-32]+'matches_found.csv')
