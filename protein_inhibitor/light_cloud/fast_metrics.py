#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
# molecular mechanics
from rdkit import Chem, RDLogger


# In[7]:


config = {"source_data" : "data/dataset.smi",
          "genereated"  : "product/generated_new.smi",
          "destin_valid": "product/valid.smi"
         }
RDLogger.DisableLog('rdApp.*')


# In[9]:


def valid(seqs):
    ''' Returns the valid smiles. '''
    valid = []
    for i, seq in enumerate(seqs): 
        mol = Chem.MolFromSmiles(seq)
        if mol is not None:
            valid.append(Chem.MolToSmiles(mol))
    return valid


def unique(seqs):
    ''' Returns the unqiue smiles. '''
    return set([seq.lower() for seq in seqs])


def original(source, seqs):
    ''' Returns the original smiles not present in the source dataset. '''
    source = [smiles.lower() for smiles in source]
    return [seq for seq in seqs if seq.lower() not in source]


def read_smiles(path):
    ''' Reads all smiles in a file. '''
    with open(config["genereated"], "r") as f:
        smiles = [line.replace("\n", "") for line in f.readlines()]
    return smiles 


def save_smiles(seqs, destin):
    ''' Saves smiles in a file. '''
    with open(config["destin_valid"], "w") as f:
        for drug in seqs: 
            f.write(drug+"\n")
        f.write("\n")


# In[12]:

source    = read_smiles(config["source_data"])
generated = read_smiles(config["genereated"])

# In[ ]:


valid = valid(generated)
# save valid smiles:
save_smiles(valid, config["destin_valid"])
print("saved valid smiles!")
# calculate validity    metric
print("Validity: {0}%".format(len(valid)/len(generated) * 100))
# calculate uniqueness  metric
print("Uniqueness: {0}%".format(len(unique(generated))/len(generated) * 100))
# calculate originality metric
print("Originality: {0}%".format(len(original(source, generated))/len(generated) * 100))

