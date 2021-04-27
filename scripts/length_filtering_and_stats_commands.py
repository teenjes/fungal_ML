
"""
The goal of this program is to examine the distribution of reads within
each file, and for the files generated after homology analysis.
The program will generate summary statistics for the result of the homology
analysis and save figures illustrating the read distribution for the
frDNA reads
"""

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord as SR
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import seaborn as sns
import matplotlib.pyplot as plt
import argparse


WHITE='\033[1;37m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'
PURPLE='\033[1;35m'

parser = argparse.ArgumentParser(description="""
The goal of this program is to examine the distribution of reads within
each file, and for the files generated after homology analysis.
The program will generate summary statistics for the result of the homology
analysis and save figures illustrating the read distribution for the
frDNA reads
""")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
parser.add_argument("full_file", help="The full, unfiltered file containing all reads for this barcode")
parser.add_argument("input_folder", help="The destination folder within which the .paf files are generated")
parser.add_argument("output_folder", help="The destination folder for any outputs from this script - including summary statistics file and plots")
args = parser.parse_args()

print('\033[0;35m'+'START'+'\033[1;37m')

output_folder = args.output_folder.rsplit('/', 1)[-2]
input_folder = args.input_folder.rsplit('/', 1)[-2]
if args.verbose:
    print('\033[0;31m' + "Input folder is " + input_folder + '\033[1;37m')
    print('\033[0;31m' + "Output folder is " + output_folder + '\033[1;37m')
    print('\033[0;34m' + "Loading " + args.full_file + '\033[1;37m')

# Load the full file containing all reads for this barcode
full_file_dict = SeqIO.to_dict(SeqIO.parse(args.full_file, "fastq"))

if args.verbose:
    print('\033[0;34m' + "Loaded " + args.full_file + '\033[1;37m')

# Extract the information about the lengths of the sequence for each read in this barcode
full_lengths = []
for key in full_file_dict:
    full_lengths.append(len(full_file_dict[key].seq))
full_lengths_len = len(full_file_dict)


# Plot the spread of read lengths for this barcode
    # Expect to see two peaks - one for EF1a and one for frDNA
ax = sns.distplot(full_lengths, color="k", kde=False, bins=5000)
ax.set(xlim=(250, 3500))
ax.set_title("Read spread for %s" % '/'.join(args.full_file.rsplit('/')[-3:-1]), fontsize=15)
ax.set_xlabel("Length of read", fontsize=13)
ax.set_ylabel("Number of reads", fontsize=13)
figure1 = ax.get_figure()
# Save this figure out
figure1.savefig('/'.join([output_folder, 'full_read_spread.png']))
figure1.clf()
if args.verbose:
    print('\033[0;32m' + "Full spread image file saved to " + '/'.join([output_folder, 'full_read_spread.png']) + '\033[1;37m')
    print('\033[0;34m' + "Loading " + input_folder+"/combined_test.paf" + '\033[1;37m')

    
    
    
    
# Import the PAF file resulting from the minimap2 homology filtering
full_paf = pd.read_csv(input_folder+"/combined_test.paf", sep='\t', header=None, engine='python')
if args.verbose:
    print('\033[0;34m' + "Loaded " + input_folder+"/combined_test.paf" + '\033[1;37m')
# Determine all the read ids present within the homology-filtered dataset
# Then, create a dictionary extracting all the information from the full read file, but ONLY for reads present within the homology-filtered data
full_dict = {}
for key in full_paf[0].unique():
    full_dict[key] = full_file_dict[key]

    
    
    
    
    
# For each key in the homology-filtered dictionary, extract the sequence length and key
full_paf_lengths = []
full_keys = []
for key in full_dict:
    full_paf_lengths.append(len(full_dict[key].seq))
    full_keys.append(key)

mean = np.mean(full_paf_lengths)
std = np.std(full_paf_lengths)

if args.verbose:
    print('\033[1;33m' + 'Mean read length is %s' % mean + '\033[1;37m')
    print('\033[1;33m' + 'Standard deviation of read length is %s' % std + '\033[1;37m')
    
    
    
    
length_filt_dict = full_dict.copy()
for key in full_keys:
    if len(full_dict[key].seq) < (mean-1.645*std) or len(full_dict[key].seq) > (mean+1.645*std):
        del length_filt_dict[key]

        
        
SeqIO.write(length_filt_dict.values(), '/'.join([output_folder, 'length_restricted_reads.fasta']), "fasta")
if args.verbose:
    print('\033[1;36m' + 'Saved %s' % ('/'.join([output_folder, 'length_restricted_reads.fasta'])) + '\033[1;37m')     
    
    
    
length_filt_lens = []
len_filt_keys = []
for key in length_filt_dict:
    length_filt_lens.append(len(length_filt_dict[key].seq))
    len_filt_keys.append(key)


    
# Create a dictionary containing the statistics for the filtered dataset
    # Total no. frDNA reads, Min. read length, Max. read length, Mean read length, Median read length

stats_dict = {'number of frDNA reads':len(length_filt_lens),'minimum read length':min(length_filt_lens),'maximum read length':max(length_filt_lens),'mean read length':"{:.0f}".format(np.mean(length_filt_lens)),'std dev':"{:.0f}".format(np.std(length_filt_lens)),'median read length':"{:.0f}".format(np.median(length_filt_lens))}

stats = pd.DataFrame(stats_dict, index=['%s' % '/'.join(args.full_file.rsplit('/')[-3:-1])])    
              
bx = sns.distplot(length_filt_lens, color="k", kde=False)
bx.set(xlim=(250, 3500))
bx.set_title("frDNA reads for %s" % '/'.join(args.full_file.rsplit('/')[-3:-1]), fontsize=15)
bx.set_xlabel("Length of read", fontsize=13)
bx.set_ylabel("Number of reads", fontsize=13)
figure2 = bx.get_figure()
figure2.savefig('/'.join([output_folder, 'frDNA_len_filt_full.png']))
figure2.clf()
if args.verbose:
    print('\033[0;32m' + "frDNA spread image file saved to " + '/'.join([output_folder, 'frDNA_len_filt_full.png']) + '\033[1;37m')

cx = sns.distplot(length_filt_lens, color="k", kde=False)
cx.set(xlim=((mean-1.645*std)-100, (mean+1.645*std)+100))
cx.set_title("frDNA reads for %s" % '/'.join(args.full_file.rsplit('/')[-3:-1]), fontsize=15)
cx.set_xlabel("Length of read", fontsize=13)
cx.set_ylabel("Number of reads", fontsize=13)
figure3 = cx.get_figure()
figure3.savefig('/'.join([output_folder, 'frDNA_len_filt_limited.png']))
figure3.clf()
if args.verbose:
    print('\033[0;32m' + "Zoomed-in frDNA spread image file saved to " + '/'.join([output_folder, 'frDNA_len_filt_limited.png']) + '\033[1;37m')

stats.to_csv('/'.join([output_folder, 'frDNA_len_filt_statistics.csv']), index=False)
if args.verbose:
    print('\033[0;32m' + "Summary statistics file saved to " + '/'.join([output_folder, 'frDNA_len_filt_statistics.csv']) + '\033[1;37m')
    
print('\033[0;35m'+'END'+'\033[1;37m')
