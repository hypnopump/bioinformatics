import numpy as np
import tensorflow as tf 

from tqdm import tqdm
from chem_utils import *
# following import only works with conda+rdkit installed:
# not essential for prediction, but for validation of smiles
try:
    from rdkit import Chem, RDLogger
except: 
    pass

class SmilesPredictor(object):
    def __init__(self, model, max_length, sampling_temp,
                 end_char, batch_predict, verbose=False):
        ''' Instantiate smiles predictor. '''
        self.st            = SmilesTokenizer()
        self.model         = model
        self.max_length    = max_length
        self.sampling_temp = sampling_temp
        self.end_char      = end_char
        self.batch_predict = batch_predict
        self.verbose       = verbose

    def _generate(self, seqs):
        ''' generates the smiles sequence. '''
        stage  = []
        indexs = list(range(len(seqs)))
        tokens = [self.st.tokenize_fast(seq) for seq in seqs]
        # predict next sample for tokens
        while len(indexs) > 0: 
            # tokenize
            tokenized = [self.st.one_hot_encode(tokens[index])
                         for i, index in enumerate(indexs)]
            # parallel prediction
            preds     = self.model.predict(np.array(tokenized))
            # sample
            next_idx  = [self.sample_with_temp(pred[-1])
                         for i,pred in enumerate(preds)]
            # add new item and check if pred is finished
            for i, index in enumerate(indexs):
                tokens[index].append(self.st.table[next_idx[i]])
                # check if smiles is done
                if (tokens[index][-1] == self.end_char) or \
                   (len(tokens[index]) > self.max_length+2):
                    stage.append(index)
            # stop elongating staged smiles
            for index in stage: 
                if index in indexs: 
                    indexs.remove(index)

        # final formatting 
        return [''.join(seq[1:]).rstrip(self.end_char) for seq in tokens]


    def sample_with_temp(self, preds):
        ''' probabilistic sampling from prediction. '''
        streched = np.log(preds) / self.sampling_temp
        # softmax function
        streched_probs = np.exp(streched) / np.sum(np.exp(streched))
        return np.random.choice(range(len(streched)), p=streched_probs)


    def check_valid(self, seqs):
        ''' Filters non-valid smiles. '''
        RDLogger.DisableLog('rdApp.*')

        checked = []
        for sample in seqs:
            mol = Chem.MolFromSmiles(seq)
            if mol is not None:
                checked.append(Chem.MolToSmiles(mol))

        return checked


    def smiles_predict(self, num=1, start='G', valid=False):
        ''' predicts "num" smiles in batches. '''
        primer = [start for i in range(self.batch_predict)]
        
        if self.verbose: 
            print("Generating smiles. Batches of {0}".format(self.batch_predict))
            sampled = [self._generate(primer) for i in\
                       tqdm(range(num//self.batch_predict + 1))]
        else:
            sampled = [self._generate(primer) for i in range(num)]
        # flatten
        sampled = [smiles for batch in sampled for smiles in batch]
        # check if valid
        return self.check_valid(sampled) if valid else sampled