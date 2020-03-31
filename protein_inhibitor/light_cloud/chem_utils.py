import copy
import json
import numpy as np
from tqdm import tqdm
try: 
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')
except:
    pass


def read_json(path):
    ''' returns a dict out of a json file. '''
    with open(path, "r") as f: 
        config = json.load(f)
    return config

def read_smiles(path):
    ''' Reads all smiles in a file. '''
    with open(path, "r") as f:
        smiles = [line.replace("\n", "") for line in f.readlines()]
    return [smile for smile in smiles if smile not in [""]]


def save_smiles(smiles, path):
    ''' Saves smiles in a file. '''
    with open(path, "w") as f:
        for smile in smiles: 
            f.write(smile+"\n")
        f.write("\n")


def save_sdf(smiles, path, mol_names=None):
    """ Writes a set of smiles in sdf format. """
    w = Chem.SDWriter(path)
    for i,smile in enumerate(smiles): 
        mol = Chem.MolFromSmiles(smile)
        mol.SetProp("_Name", str(mol_names[i]))
        w.write(mol)


def save_sdf_ind(smiles, path, mol_names=None):
    """ Writes a N smiles in N separate PDB files. """
    for i,smile in enumerate(smiles):
        w = Chem.SDWriter(path+str(mol_names[i])+".sdf") 
        mol = Chem.MolFromSmiles(smile)
        mol.SetProp("_Name", str(mol_names[i]))
        w.write(mol)


def save_pdb(smiles, path, mol_names=None):
    """ Writes a N smiles in N separate PDB files. """
    for i,smile in enumerate(smiles):
        w = Chem.PDBWriter(path+str(mol_names[i])+".pdb") 
        mol = Chem.MolFromSmiles(smile)
        mol.SetProp("_Name", str(mol_names[i]))
        w.write(mol)


def validity(seqs):
    ''' Returns the valid smiles. '''
    valid = []
    for seq in seqs: 
        mol = Chem.MolFromSmiles(seq)
        if mol is not None:
            valid.append(Chem.MolToSmiles(mol))
    return valid


def uniqueness(seqs):
    ''' Returns the unqiue smiles. '''
    return set([seq for seq in seqs])


def originality(source, seqs, process_source=None, verbose=None):
    ''' Returns the original smiles not present in the source dataset. '''
    if process_source:
        source = [Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
                  for smiles in source]
    if verbose:
        return [seq for seq in tqdm(seqs) if seq not in source]
    else: 
        return [seq for seq in seqs if seq not in source]

        
class SmilesTokenizer(object):
    def __init__(self):
        atoms = [
            'Li',
            'Na',
            'Al',
            'Si',
            'Cl',
            'Sc',
            'Zn',
            'As',
            'Se',
            'Br',
            'Sn',
            'Te',
            'Cu',
            'Mg',
            'Mn',
            'Ca',
            'Fe',
            'Te',
            'Co',
            'H',
            'B',
            'C',
            'N',
            'O',
            'F',
            'P',
            'S',
            'K',
            'V',
            'I',
        ]
        special = [
            '(', ')', '[', ']', '=', '#', '%', '0', '1', '2', '3', '4', '5',
            '6', '7', '8', '9', '+', '-', 'se', 'te', 'c', 'n', 'o', 's', '/',
            '\\', '@', '@@', '.', 'i', 'p', 'f', 's' 
        ]
        padding = ['G', 'A', 'E']

        self.table = sorted(atoms, key=len, reverse=True) + special + padding
        self.table_len = len(self.table)

        self.one_hot_dict = {}
        for i, symbol in enumerate(self.table):
            vec = np.zeros(self.table_len, dtype=np.float32)
            vec[i] = 1
            self.one_hot_dict[symbol] = vec

    def tokenize(self, smiles):
        N = len(smiles)
        i = 0
        token = []
        while (i < N):
            for j, symbol in enumerate(table):
                if symbol == smiles[i:i + len(symbol)]:
                    token.append(symbol)
                    i += len(symbol)
                    break
        return np.array(token)

    def tokenize_fast(self, smiles):
        N = len(smiles)
        i = 0
        token = []
        while (i < N):
            token_a, token_b = smiles[i:i+2], smiles[i:i+1]
            # first try 2 characters, then 1
            if token_a in self.table:
                token.append(token_a)
                i += 2
            elif token_b in self.table:
                token.append(token_b)
                i += 1
            else:
                print("token not found @ pos:", i, "in", smiles[:100])
                i +=1

        return token

    def one_hot_encode(self, tokenized_smiles):
        result = np.array(
                    [self.one_hot_dict[symbol] for symbol in tokenized_smiles],
                 dtype=np.float32)
        # result = result.reshape(1, result.shape[0], result.shape[1])
        return result