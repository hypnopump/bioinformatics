""" Trying to determine whether a molecule is a 
	valid smiles encoding or not. 
	Standardize smiles in order to check if 2 encodings
	are the same molecule.
"""
import numpy as np 
from chempy import Substance


vocab = ["H", "Li", "Be", "B", "C", "N", "O", "F", "=", "#", "(", ")", "1", "2", "3"]
vocab_sec = ["Na", "Mg", "Al", "Si", "P", "S", "Cl", "K", "Ca", "Fe" "Cu", "Zn", "Se", "Te", "Br", "I"]

a = "CN1C=NC2=C1#C(#O)N(C(=O)N2CC"
b = Substance.from_formula(a)

print(a, b)



