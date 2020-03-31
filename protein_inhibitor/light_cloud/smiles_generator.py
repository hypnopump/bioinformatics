import numpy as np
import tensorflow as tf
# Custom functions import
from chem_utils import *


class SmilesGenerator(tf.keras.utils.Sequence):
    ''' Generates data for the NN '''
    def __init__(self, smiles, batch_size=128, shuffle=True,
                       max_length=256, verbose=False):
        ''' Initialization. '''
        self.st         = SmilesTokenizer()
        self.max_length = max_length
        self.batch_size = batch_size
        self.shuffle    = shuffle
        self.smiles     = smiles
        self.indexes    = np.arange(len(self.smiles))
        self.on_epoch_end()


    def __len__(self):
        ''' Denotes the number of batches per epoch. '''
        return int(np.floor(len(self.smiles) / self.batch_size))


    def __getitem__(self, index):
        ''' Generate one batch of data. '''
        indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]
        x, y    = self.__data_generation(indexes)

        return x, y


    def on_epoch_end(self):
        ''' Updates indexes after each epoch. '''
        if self.shuffle == True:
            np.random.shuffle(self.indexes)
                

    def __data_generation(self, indexes):
        ''' Generates data containing batch_size samples. '''
        tokenized = [self.st.tokenize_fast(self.smiles[i]) for i in indexes]
        padded    = [["G"]+smile+["A"]+["E" for i in range(self.max_length-len(smile)-2)]
                     for smile in tokenized]
        onehot    = np.array([self.st.one_hot_encode(smile) for smile in padded])

        return onehot[:, :-1, :], onehot[:, 1:, :]