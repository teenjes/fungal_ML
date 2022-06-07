
"""
Script aims to read in the reference_dataframe file, select a taxonomic
level and group, and read the path to the location of that data. It then
prepares data for machine learning by converting base pair coding to numerical
encoding, pads it out and then runs the algorithm
"""

## REQUIRED PACKAGES
# argparse, Biopython, keras, math, matplotlib, numpy, os, pandas, random, scikit-learn

import pandas as pd
from Bio import SeqIO
import numpy as np
import os
import random
import argparse
import keras
from keras.models import Sequential
from keras.layers import Dense, Conv1D, Dropout, MaxPooling1D, Flatten
from keras.utils.vis_utils import plot_model
import math
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix,classification_report,accuracy_score,precision_score,recall_score,f1_score


## Definitions of functions

# Return the longest sequence length of filtered sequences
def max_seq_len(SeqIO_dict):
    """
    Function takes a SeqIO_dict and returns the lengths of the
    longest sequence
    """
    total_lens = []
    for key in SeqIO_dict.keys():
        total_lens.append(len(SeqIO_dict[key].seq))
    return max(total_lens)

# Convert nucleic data to numerical data and pad sequence to pad_length
def numberfy(SeqIO_dict, seq_len, nsubsample, genus_name, species_name, pad_length):
    """
    Take SeqIO_dict and return SeqIO_dict were bases have been replaced
    with numbers
    ACGT- replaced with 01234
    Take the seq_len each sequence should have
    """
    num_dict = {}
    
    keys = list(SeqIO_dict.keys())
    randkeys = random.sample(keys, k=nsubsample)
    
    with open(data_root+'models/ids_%s_%s_%s_%s_%s.txt' % (args.tax_rank,args.name,args.n_reads,genus_name,species_name), 'w+') as file:
        file.writelines("%s\n" % key for key in randkeys)
    
    for key in randkeys:
        seq = str(SeqIO_dict[key].seq).replace("A",'0 ')\
        .replace("C",'1 ').replace("G",'2 ').replace("T",'3 ')\
        .replace("a",'0 ').replace("c",'1 ').replace("g",'2 ')\
        .replace("t",'3 ')
#         seq_new = seq + '4 '*(seq_len - int(len(seq)/2))
        seq_new = seq + '4 '*(pad_length - int(len(seq)/2))
        if seq_new.find('t') != -1:
            print(seq_new.find('t'))
            print("ERROR - strange value in sequence")
            print(seq_new)
            exit()
        num_dict[key] = list(map(int, seq_new.split(' ')[:-1]))
    return num_dict

## Define the variables to allow the script to be run on the command line with variable input
parser = argparse.ArgumentParser(description="""
Script aims to read in the reference_dataframe file, select a taxonomic
level and group, and read the path to the location of that data. It then
prepares data for machine learning by converting base pair coding to numerical
encoding, pads it out and then runs the algorithm
""")
parser.add_argument("ref_df_fn", help="File path to the reference dataframe")
parser.add_argument("data_root", help="Root folder for analysis/")
parser.add_argument("--tax_rank", "-r", help="taxonomic rank for analysis (e.g. phylum, genus)")
parser.add_argument("--name", "-n", help="name of rank to select from (e.g. ascomycota, candida)")
parser.add_argument("--n_reads", "-c", help="count of reads per class (e.g. 15000)")
parser.add_argument("--len", "-l", help="length to pad sequences out to (e.g. 5000)")
group = parser.add_mutually_exclusive_group()
group.add_argument("--verbose", "-v", "--v", action="store_true")
group.add_argument("--quiet", "-q", "--q", action="store_true")
args = parser.parse_args()

# assign required arguments to variables
ref_df_fn = args.ref_df_fn
data_root = args.data_root

# assign a number of reads per class
n_reads = int(args.n_reads)

# Define a padded length if one is not specified
if args.len:
    pad_len = args.len
else:
    pad_len = 5000

# test to make sure both required file paths are input
try:
    os.path.exists(ref_df_fn)
except:
    print('Cannot find %s' % ref_df_fn)
try:
    os.path.exists(data_root)
except:
    print('Cannot find %s' % data_root)
    
if args.verbose:
    print('\033[1;34m' + "Reference dataframe is at " + ref_df_fn + '\033[0m')
    print('\033[1;34m' + "Root directory is at " + data_root + '\033[0m')
    if args.tax_rank and args.name:
        print('\033[1;34m' + "Tax Rank is " + args.tax_rank.lower() + '\033[0m')
        print('\033[1;34m' + "Name is " + args.name.lower() + '\033[0m')
    elif args.one and args.two:
        print('\033[1;34m' + "Species one is " + args.one.lower() + '\033[0m')
        print('\033[1;34m' + "Species two is " + args.two.lower() + '\033[0m')
    print('\033[1;34m' + "Count of reads per sample is", n_reads,'\033[0m')
    
    
    
# read in the reference dataframe from the argument path
ref_df = pd.read_csv(ref_df_fn, index_col=None)

# check whether the reference dataframe implies there are enough reads
# to continue given n_reads
try:
    if ref_df[ref_df["# for use"] \
              < n_reads].shape[0] > 0 :
        print("These species need more reads.")
        print(ref_df[ref_df["# for use"] \
              < n_reads])
        exit()
except:
    print('Check %s to have the wanted column names' % ref_df_fn)
    
# assign flagged variables as lower case and assign indices dataframe
if args.tax_rank and args.name:
    tax_rank = args.tax_rank.lower()
    name = args.name.lower()
    try:
        indices = ref_df[ref_df[tax_rank] == name].index
        print(indices)
    except:
        print("Tax_rank or Name not found in reference_dataframe")
        print(tax_rank, name)
elif args.one and args.two:
    one = args.one.lower()
    two = args.two.lower()
    try:
        indices = ref_df[(ref_df['species'] == one)&(ref_df['species'] == two)].index
        print(indices)
    except:
        print("Species inputs not found in reference_dataframe")
        print(one, two)

# where the values are that index's path's dataframe
SeqIO_dicts = {}
for index in indices:
    fasta_path = ref_df.loc[index, 'path for use']
    try:
        SeqIO_dicts[index] = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    except:
        print('Check location of fasta files')
        print(fasta_path, "does not exist")
        
# each path within an index corresponds to a species
# if tax_rank > genus, we want to look at which species are within which genus/family/order etc.

# determine the maximum sequence length of accepted sequences
total_lens = []
for key, value in SeqIO_dicts.items():
    total_lens.append(max_seq_len(value))
print('\033[0;32m'+"The maximum sequence length of all sampled sequences is"+ '\033[1;37m',max(total_lens),'\033[0m')


# randomly subsample n_reads number of reads from each index's corresponding
# set of reads, convert base pair coding to numerical coding and 
# pad to the max sequence length
numSeqIO_dicts = {}
max_len = max(total_lens)
del total_lens
if (args.one and args.two) or tax_rank == "genus":
    for key, value in SeqIO_dicts.items():
        numSeqIO_dicts[key] = numberfy(value, max_len, n_reads,ref_df.loc[key,'genus'], ref_df.loc[key,'species'], pad_len)
else:
    location = (ref_df.columns.get_loc(tax_rank)-1)
    col_name = ref_df.columns[location]
    if args.verbose:
        print('location is', col_name)

    classes = ref_df.iloc[indices,location].unique()
    if args.verbose:
        print('classes are', classes)

    count_dict = {}
    for class_ in classes:
        count_dict[class_] = sum(ref_df.iloc[indices,location] == class_)
    if args.verbose:
        print('count_dict is', count_dict)
    
    del classes
        
    min_vals = []
    for class_, n_class in count_dict.items():
#         if n_class == min(count_dict.values()):
            min_vals.append(ref_df[ref_df.iloc[:,location] == class_]['# for use'].min())
    if min(min_vals) % 2 == 0:
        minimum_value = int(min(min_vals))
    else:
        minimum_value = int(min(min_vals)-1)
    if minimum_value > 35000:
        minimum_value = 35000
        
    del min_vals
    
    if args.verbose:
        print('minimum number of reads is', minimum_value)
    class_lens_ind = []
    if len(count_dict) > 1:
        max_reads = 0
        for key, value in count_dict.items():
            if value == max(count_dict.values()):
                max_reads = value*n_reads

        if max_reads <= minimum_value:
            minimum_value = max_reads

        for key, n_class in count_dict.items():
            s_reads = int(minimum_value/n_class)
            if ref_df[ref_df.loc[:,col_name]==key]['# for use'].min() < s_reads:
                minimum_value = ref_df[ref_df.loc[:,col_name]==key]['# for use'].min()/n_class
                s_reads = int(minimum_value/n_class)
            if args.verbose:
                print('The class is', key, 'and the number of reads to be subsampled is', s_reads)
            for keya, value in SeqIO_dicts.items():
                if ref_df.loc[keya,col_name] == key:
                    numSeqIO_dicts[keya] = numberfy(value, max_len, s_reads, ref_df.loc[keya,'genus'], ref_df.loc[keya,'species'], pad_len)
                    class_lens_ind.append(s_reads)
        n_reads = minimum_value
    elif len(count_dict) == 1:
        s_reads = n_reads
        print("no comparison for the rank")
        exit()

        
try:
    del count_dict
except:
    pass


location = (ref_df.columns.get_loc(tax_rank)-1)
col_name = ref_df.columns[location]
classes = ref_df.iloc[indices,location].unique()

del indices

order = []
seq_list = []
total_expected_reads = len(classes)*n_reads
class_lens = []
for class_ in classes:
    tmp_sum = []
    for key in numSeqIO_dicts.keys():
        if ref_df.loc[key,col_name] == class_:
            order.append(key)
            seq_list.append(np.array(list(numSeqIO_dicts[key].values())))
            tmp_sum.append(len(list(numSeqIO_dicts[key].values())))
    class_lens.append(sum(tmp_sum))

try:
    del tmp_sum
except:
    pass
try:
    del numSeqIO_dicts
except:
    pass
    
total_actual_reads = min(class_lens)

print(class_lens)
if args.verbose:
    print("Ids order for labels is", order)
    if not (args.one and args.two) and tax_rank != "genus":
        print("Number of reads subsampled per id is", class_lens_ind)
    print("Total expected reads is", total_expected_reads)
    for i in range(0, len(classes)):
        print(classes[i], "has", class_lens[i], "reads")
    print("Total reads used per class is", sum(class_lens))
    print("Total actual reads available per class is", total_actual_reads)

try:
    del class_lens_ind
except:
    pass
    
seq_comb = np.concatenate(seq_list, axis = 0)
num_class = len(classes)

try:
    del seq_list
except:
    pass

if len(set(class_lens)) == 1:
    all_data = seq_comb
else:
    class_lens_cumsum = np.cumsum(class_lens)
    new_seq_list = []
    for i in range(0, len(class_lens_cumsum)):
        if i == 0:
            new_seq_list.append(seq_comb[0:class_lens_cumsum[i]][:total_actual_reads])
        else:
            new_seq_list.append(seq_comb[class_lens_cumsum[i-1]:class_lens_cumsum[i]][:total_actual_reads])
    del seq_comb
    all_data = np.concatenate(new_seq_list, axis = 0)

try:
    del class_lens
except:
    pass
try:
    del new_seq_list
except:
    pass
    
    
    
    
# determine the number of classes and generate an array of ids
all_labels_onehot = np.zeros( (total_actual_reads*num_class,num_class) )
for i in range(0, num_class):
    all_labels_onehot[i*total_actual_reads:(i+1)*total_actual_reads,i] = 1

# Print the shape of the resulting dataframes to visually verify
if args.verbose:
    print('all_labels_onehot.shape: ', all_labels_onehot.shape)
    print('all_data.shape:', all_data.shape)

samples_count = total_actual_reads*num_class
if args.verbose:
    print('samples_per_class:', total_actual_reads)
    print('samples_count:', samples_count)

# Create a method for shuffling data
shuffle_indices = random.sample(range(0, samples_count), samples_count)
if args.verbose:
    print(len(shuffle_indices))

# Assign a percentage of data for training and the rest for testing and define which reads belong to which dataset
train_size = math.floor(0.85*all_data.shape[0])
if args.verbose:
    print("Training data size:", train_size)
indices_train = shuffle_indices[0:train_size]
indices_test = shuffle_indices[train_size+1:samples_count]

try:
    del train_size
except:
    pass
try:
    del samples_count
except:
    pass
try:
    del shuffle_indices
except:
    pass

# Define the data vs labels for each of the training and test sets
X_train = all_data[indices_train,:]
Y_train = all_labels_onehot[indices_train]
X_test = all_data[indices_test,:]
Y_test = all_labels_onehot[indices_test]

try:
    del all_data
except:
    pass

if args.verbose:
    print('X_train.shape : ', X_train.shape)
    print('X_test.shape : ', X_test.shape)
    print('Y_train.shape : ', Y_train.shape)
    print('Y_test.shape : ', Y_test.shape)

# Define the input dimension from X_train.shape[1]
in_dim = X_train.shape[1]

classes_dict = {}
for i in range(0, len(classes)):
    classes_dict['%s' % i] = classes[i]

    
## Define and train the model
model = Sequential()
model.add(Conv1D(filters=64, kernel_size=3, activation='relu', input_shape=(X_train.shape[1],1)))
model.add(Conv1D(filters=64, kernel_size=3, activation='relu'))
model.add(Dropout(0.5))
model.add(MaxPooling1D(pool_size=2))
model.add(Flatten())
model.add(Dense(100, activation='relu'))
model.add(Dense(num_class, activation='softmax'))
model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
model.summary()
print()
history = model.fit(np.expand_dims(X_train,2), Y_train, validation_data=(np.expand_dims(X_test,2), Y_test), batch_size=128, epochs=10, verbose=1)


# Save the history of the model training plot for accuracy vs val_accuracy
model.save(data_root+'models/model_%s_%s_%s.h5' % (args.tax_rank,args.name,args.n_reads))
plt.plot(history.history['accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('Model accuracy for %s %s' % (name, tax_rank))
plt.ylabel('Accuracy')
plt.xlabel('Epoch')
plt.legend(['Train', 'Test'], loc='upper left')
# plt.show()
plt.savefig(data_root+'plot_histories/history_%s_%s_%s.png' % (args.tax_rank,args.name,args.n_reads))
plt.close()

# Define the predicted output for the test set
yhat_probs = model.predict(np.expand_dims(X_test,2), verbose=0)
print(yhat_probs.shape)
print(yhat_probs)
yhat_classes = np.argmax(yhat_probs,axis=1)
print(yhat_classes.shape)
print(yhat_classes)

Y_test_ints = np.where(Y_test==1)[1]
print(Y_test_ints.shape)
yhat_probs = yhat_probs[:, 0]
yhat_classes = yhat_classes[:]

## Using the test dataset, measure accuracy metrics for the model

# accuracy: (tp + tn) / (p + n)
accuracy = accuracy_score(Y_test_ints, yhat_classes)
print('Accuracy: %f' % accuracy)

# precision tp / (tp + fp)
precision = precision_score(Y_test_ints, yhat_classes, average=None)
print('precision: ', precision)

# recall: tp / (tp + fn)
recall = recall_score(Y_test_ints, yhat_classes, average=None)
print('recall: ', recall)

# f1: 2 tp / (2 tp + fp + fn)
f1 = f1_score(Y_test_ints, yhat_classes, average=None)
print('f1: ', f1)

# confusion matrix
matrix = confusion_matrix(Y_test_ints, yhat_classes)
print(matrix)
crosstab = pd.crosstab(Y_test_ints, yhat_classes, rownames=['True'], colnames=['Predicted'], margins=True)
print(crosstab)

# Save these statistics and confusion matrix to an output csv file
data = [accuracy]
datacol =['accuracy']
count = 1
for i in precision:
    data.append(i)
    datacol.append('precision%i'%count)
    count += 1
count = 1
for i in recall:
    data.append(i)
    datacol.append('recall%i'%count)
    count += 1
count = 1
for i in f1:
    data.append(i)
    datacol.append('f1%i'%count)
    count += 1
stats = pd.DataFrame(data=[data],columns=datacol)
stats.to_csv(data_root+'models/stats_%s_%s_%s.csv' % (args.tax_rank,args.name,args.n_reads))
crosstab.to_csv(data_root+'models/confusion_%s_%s_%s.csv' % (args.tax_rank,args.name,args.n_reads))

# Save a csv that details how the output classes of the model (e.g. 0, 1) correspond to the taxonomic output (e.g. ascomycota, basidiomycota)
with open(data_root+'models/keys_%s_%s_%s.csv' % (args.tax_rank,args.name,args.n_reads), 'w+') as f:
    for key in classes_dict.keys():
        f.write("%s,%s\n"%(key,classes_dict[key]))
if args.verbose:
    print(classes_dict)
