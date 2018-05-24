# GAN to determine validity of smiles

import numpy as np 
import pandas as pd 

origin = "data.txt"
vocab_list = [' ', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'n', '@', '-', 'AL', ')', '[', ']', 'Na', 'Ba', 'Be',
		 	  'Mg', 'Ca', 'K', 'Br', 'C', 'Te', 'I', 'S', 'Se', 'Si', 's', 'se', 'si', 'o', 'b', 'br', 'Li', 'c', 'B',
		 	  '.', '=', '+', 'N', 'F', 'Fe', 'Cl', '\\', '/', '\t', '#', '(', 'l', '5', 'P', 'H', 'O']

vocab_special = [c for c in vocab_list if len(c)>1]
print(vocab_special)

# Create vocab encoder and decoder
vocab_enc = {}
for i, el in enumerate(vocab_list):
	vocab_enc[el] = i
vocab_dec = {}
for i, el in enumerate(vocab_list):
	vocab_dec[i] = el

def extract(filename):
	""" Extract data from the txt file."""
	t = pd.read_csv('data.txt', header = None).values.tolist()
	return [subs[0] for subs in t]

data = extract(origin)
print(data[:5])

print(vocab_list)



n_compounds = len(data)
n_vocab = len(vocab_list)
print("Total Compounds for training: ", n_compounds)
print("Total Vocab: ", n_vocab)

def encode(compound):
	code = []
	while len(compound)>0:
		if compound[:2] in vocab_special:
			code.append(vocab_enc[compound[:2]])
			compound = compound[2:]
		else:
			code.append(vocab_enc[compound[:1]])
			compound = compound[1:]
	return code

# print(encode(data[0]))

def data_prep(data, seq_len=60): 
	# measure max length - 60
	# seq_len = max([len(encode(subs)) for subs in data])
	data = [encode(subs) for subs in data]
	datex = []
	for subs in data:
		while len(subs)<60:
			subs.append(0)
		datex.append(subs)
	return np.array(datex)

data_ready = data_prep(data)
print("READY DATA")
print(data_ready[:3])
print(data_ready.shape)

