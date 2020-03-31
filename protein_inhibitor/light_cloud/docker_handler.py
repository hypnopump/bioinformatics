#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import warnings
warnings.filterwarnings("ignore")

# science
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# deep learning
import tensorflow as tf
# molecular mechanics
try:
    from rdkit import Chem, RDLogger
except:
    pass
# utils
from tqdm import tqdm
from chem_utils import *
from smiles_generator import *
from smiles_predictor import *


class DockerHandler():
	def __init__(self, openbabel_path="obabel", vina_path="vina",
				 obminimize_path="obminimize", docking_config=False,
				 docking_receptor=False):
		''' Handles the docking and energy scoring for ligand+receptor.
			Uses:
			* Obminimize (energy minimization)
			* OpenBabel (.pdbqt conversion)
			* AutoDock Vina (molecular docking) 
			* Optional:
				* docking_config: str. config path for vina.
				* docking_receptor: sr. receptor path for vina.
		'''
		self.openbabel_path  = openbabel_path
		self.obminimize_path = obminimize_path
		self.vina_path       = vina_path
		# additional params
		self.docking_config  = docking_config
		self.docking_receptor= docking_receptor


	def minimize_energy_mol(self, path):
		''' Minimizes the energy for a molecule using Obminimize. 
			Inputs: 
			* path: str. molecule path.
		'''
		os.system('"{obminimize}" -c 1e-4 -ff Ghemical {path}'.format(
				   obminimize=self.obminimize_path, path=path))


	def pdbqt_convert(self, in_file, out_file):
		''' Converts to .pdbqt file using OpenBabel 
			Inputs: 
			* in_file: str. input molecule path.
			* out_file: str. output molecule path.
		'''
		os.system('"{babel}" {in_file} -O {out_file} --gen3d'.format(
				   babel=self.openbabel_path, in_file=in_file,
				   out_file=out_file))


	def molecular_docking(self, config, receptor, ligand, out_file, log_file, cpu=False):
		''' Converts to .pdbqt file using OpenBabel 
			Inputs: 
			* config: str. config path (.txt)
			* receptor: str. receptor path (.pdbqt)
			* ligand: str. ligand path (.pdbqt)
			* out_file: str. output path. (.pdbqt)
			* log_file: str. logs path.
			* cpu: int. Number of cpus to use. 
		'''
		cmd = '"{vina}" --config {config} --out {out_file} --log {log_file} --receptor {receptor}  --ligand {ligand}'.format(
				vina=self.vina_path, config=config, out_file=out_file,
				log_file=log_file, receptor=receptor, ligand=ligand)
		if cpu: 
			cmd = cmd+" --cpu {cpu}".format(cpu=cpu)
		os.system(cmd)
 

	def parse_vina_results(self, paths):
		''' Parses the output of Vina from multiple dokings. 
			Inputs: 
			* paths: list of vina output files to parse
			Outputs:
			* dict containing pairs of {name: energy}
		'''
		parser_dict = {"name": [], "energy":[]}
		for path in paths:
			# handle non-existence
			try: 
				with open(path, "r") as f:
					lines = f.readlines()
			except:
				continue
			# select only important rows
			lines    = [line[:-2].split(" ") for line in lines]
			selected = []
			for line in lines:
				selected.append([x for x in line if x != ""])
			for sel in selected:
				# add to parser
				try:
					n, aff      =  int(sel[0]), float(sel[1])
					ligand_name = path.split("/")[-1].split("_")[0]
					parser_dict["name"].append(int(ligand_name))
					parser_dict["energy"].append(aff)
				except:
					pass
		return parser_dict 