import numpy as np
from keras.models import Sequential
from keras.layers import Dense
from PIL import Image
import random
import math

data_npz = np.load('../analysis/arrays_test/20171103_FAH15473_b2+b6_seqs.csv.npz')
data = data_npz['arr_0']

labels_npz = np.load('../analysis/arrays_test/20171103_FAH15473_b2+b6_ids.csv.npz')
labels = labels_npz['arr_0']

print(labels.shape)
print(data.shape)
data_class_2 = data[(labels == 2)]
data_class_6 = data[(labels == 6)]
print(data_class_2[50])
print(data_class_6[50])
print('data_class_2.shape: ', data_class_2.shape)
print('data_class_6.shape: ', data_class_6.shape)
samples_per_class = data_class_2.shape[0]
samples_count = samples_per_class*2
print('samples_per_class: ', samples_per_class)
print('samples_count: ', samples_count)
all_data = np.vstack((data_class_2, data_class_6))
print('all_data.shape : ', all_data.shape)
all_labels = np.hstack( (np.zeros(samples_per_class), np.ones(samples_per_class)) )
print('all_labels.shape : ', all_labels.shape)
shuffle_indices = random.sample(range(0, samples_count), samples_count)
print(len(shuffle_indices))
train_size = math.floor(0.85*all_data.shape[0])
indices_train = shuffle_indices[0:train_size]
indices_test = shuffle_indices[train_size+1:samples_count]
X_train = all_data[indices_train,:]
Y_train = all_labels[indices_train]
X_test = all_data[indices_test,:]
Y_test = all_labels[indices_test]
print('X_train.shape : ', X_train.shape)
print('X_test.shape : ', X_test.shape)
print('Y_train.shape : ', Y_train.shape)
print('Y_test.shape : ', Y_test.shape)
# define the keras model
#model = Sequential()
#model.add(Dense(128, input_dim=480, activation='relu'))
#model.add(Dense(64, activation='relu'))
#model.add(Dense(32, activation='relu'))
#model.add(Dense(2, activation='softmax'))
model = Sequential()
model.add(Dense(32, activation='relu', input_dim=3807))
model.add(Dense(16, activation='relu'))
model.add(Dense(1, activation='sigmoid'))
model.compile(optimizer='rmsprop', loss='binary_crossentropy', metrics=['accuracy'])
# compile the keras model
#model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
model.fit(X_train, Y_train, validation_data=(X_test, Y_test), batch_size=100, epochs=40, verbose=1)
