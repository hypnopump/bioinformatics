#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import time
# science
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# deep learning
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
# molecular mechanics
try:
	from rdkit import Chem, RDLogger ,DataStructs
	from rdkit.Chem import Descriptors
except:
    pass
# big data & parallel workloads
from tqdm import tqdm
from joblib import Parallel, delayed
# custom utils
from chem_utils import *
from smiles_generator import *
from smiles_predictor import *
from model_handler import *
from docker_handler import *

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


class GensOptimizer():
	def __init__(self, config, model_handler, docker_handler, data_source,
				 base_sim=0.7, first_smiles=10000, gen_smiles=500, retain=50,
				 start_smiles=False, start_gen=0, n_start_gen=0,
				 parallel=False, n_cpu=1, verbose=False):
		''' Iteratively optimizes the produced smiles to incresae affinity
			for the desired target.
		'''
		self.config        = config
		self.model_handler = model_handler
		self.docker_handler= docker_handler
		self.data_source   = read_smiles(data_source)
		self.base_sim      = base_sim
		self.first_smiles  = first_smiles
		self.gen_smiles    = gen_smiles
		self.retain        = retain
		self.start_smiles  = start_smiles
		self.start_gen     = start_gen
		self.n_start_gen   = n_start_gen
		self.parallel      = parallel
		self.n_cpu         = n_cpu
		self.verbose       = verbose
		# additional vars
		self.params_gen    = {"batch_size": config["batch_size"], 
							  "max_length": config["smiles_max_length"],
							  "shuffle"   : config["datagen_shuffle"],
							  "verbose"   : config["datagen_verbose"] }
		# infrastructure
		self.historic_best = []
		if isinstance(self.start_gen, str):
			self.gen_counter   = self.n_start_gen
			self.index_table   = pd.read_csv(self.start_gen, sep=",")
			self.skip_eval     = True
		else: 
			self.gen_counter   = 0
			self.index_table   = None
			self.gen           = self.first_generate()
			self.skip_eval     = False


	def purify_smiles(self, smiles):
		''' Calculates stats and returns valid, unique and
			original smiles only.
		'''
		valid    = validity(smiles)
		unique   = uniqueness(valid)
		if self.index_table is None:
			source = self.data_source
		else: 
			source = self.data_source+list(self.index_table["smiles"].values)
		original = originality(source, unique)
		print(time.time(), "Starting generation #{0}".format(self.gen_counter))
		print("Begin Stats || N ={0}".format(len(smiles)))
		print("Validity = {0}%".format(100*len(valid)/len(smiles)))
		print("Uniqueness = {0}%".format(100*len(unique)/len(valid)))
		print("Originality = {0}%".format(100*len(original)/len(unique)))
		return original


	def expand_table(self, smiles, names):
		''' Adds a set of smiles to tracking table. '''
		table           = pd.DataFrame()
		table["smiles"] = smiles
		table["name"]   = names
		# affinities are expressed as kcal/mol
		table["best"]   = [None for i in range(len(smiles))]
		table["mean"]   = [None for i in range(len(smiles))]
		table["gen_cr"] = [self.gen_counter for i in range(len(smiles))]
		table["gen_ev"] = [None for i in range(len(smiles))]
		# scores
		table["w_adj"]  = [None for i in range(len(smiles))]
		table["s_adj"]  = [None for i in range(len(smiles))]
		table["adj"]    = [None for i in range(len(smiles))]
		# add to main index table or expand it
		if self.index_table is None:
			self.index_table = table
		else: 
			self.index_table = pd.concat([self.index_table, table])



	def first_generate(self):
		''' Generates the first gen of smiles. '''
		if self.start_smiles:
			smiles = read_smiles(self.start_file)
		else: 
			smiles = self.model_handler.predictor.smiles_predict(
						  num=self.first_smiles, valid=False)
		# purify
		smiles = self.purify_smiles(smiles)
		# structure smiles
		self.expand_table(smiles, np.arange(len(smiles)))
		return smiles
	

	def new_gen(self):
		''' Generates the first gen of smiles. '''
		smiles = self.model_handler.predictor.smiles_predict(
					  num=self.gen_smiles, valid=False)
		# purify
		smiles = self.purify_smiles(smiles)
		# structure smiles
		self.expand_table(smiles, np.arange(len(self.index_table),
									 len(self.index_table)+len(smiles)))


	def select_candidates(self):
		''' Selects the ideal amount to evaluate before going to next gen. '''
		gen = self.index_table[self.index_table["gen_cr"] == self.gen_counter]\
				  .values
		# newly_selected 
		smiles, names = [], []
		# handle verbosity
		rest_iterator = tqdm(gen[:, 0]) if self.verbose else gen[:, 0]
		# Prepare fingerprints for similarity calcs
		best_fingerprints    = [Chem.RDKFingerprint(Chem.MolFromSmiles(mol))
								for mol in self.historic_best]
		rest_fingerprints    = [Chem.RDKFingerprint(Chem.MolFromSmiles(mol))
		            			for mol in rest_iterator]
		# start adding non-similar structs in similarity order
		sim_threshold = 0.05 
		sim_target    = self.base_sim if self.gen_counter == 0 else 1
		while len(names) < 2500:
			for (mol, i), fingerprint in zip(gen[:, :2], rest_fingerprints):
				# find maximum similarity (1max-0min) to list of best molecs
				if len(names) == 0:
					max_sim = 0 
				else: 
					max_sim = np.max(DataStructs.BulkTanimotoSimilarity(
									 fingerprint, best_fingerprints))

				if i not in names and max_sim <= sim_threshold:
					best_fingerprints.append(fingerprint)
					smiles.append(mol)
					names.append(i)
			# logs
			print("Complete threshold: {sim_thresh} with n={n}. Target={target}"\
				  .format(sim_thresh=sim_threshold, n=len(names),
				  		  target=len(gen)))
			# break if sim_target is reached
			if sim_threshold >= sim_target:
				break 
			# increase threshold + ensure max_sim
			sim_threshold += 0.05
			if sim_threshold > 1.0: 
				sim_threshold = 1
		return smiles, names


	def evaluate_pop(self, gen, names, work_path, docked_path):
		''' Updates the index_table with new information from the population.
			Inputs:
			* work_path: working path to save .sdf and .pdbqt ligand files
			* docked_path: path to save docking logs and .pdbqt docked structs
		'''
		# save as .sdf files
		names = [str(x)+"_gen_"+str(self.gen_counter) for x in names]
		save_sdf_ind(gen, path=work_path, mol_names=names)
		# minimize
		iterator = tqdm(range(len(names))) if self.verbose else range(len(names))
		if self.parallel: 
			_ = Parallel(n_jobs=-1)\
				(delayed(self.docker_handler.minimize_energy_mol)\
						(work_path+names[i]+".sdf") for i in iterator)
		else:
			for i in iterator:
				self.docker_handler.minimize_energy_mol(work_path+\
														names[i]+".sdf")
		print(time.time(), "Energy minimization done")
		# pdbqt
		iterator = tqdm(range(len(names))) if self.verbose else range(len(names))
		if self.parallel: 
			_ = Parallel(n_jobs=-1)\
				(delayed(self.docker_handler.pdbqt_convert)\
						(work_path+names[i]+".sdf", work_path+names[i]+".pdbqt")\
						for i in iterator)
		else:
			for i in iterator:
				self.docker_handler.pdbqt_convert(work_path+names[i]+".sdf",
												  work_path+names[i]+".pdbqt")
		print(time.time(), "PDBQT conversion done")
		# dock
		iterator = tqdm(range(len(names))) if self.verbose else range(len(names))
		if False: # self.parallel: 
			_ = Parallel(n_jobs=-1)\
				(delayed(self.docker_handler.molecular_docking)\
						(self.docker_handler.docking_config,
						 self.docker_handler.docking_receptor,
						 work_path+names[i]+".pdbqt",
						 docked_path+names[i]+"_docked.pdbqt",
						 docked_path+names[i]+"_logs.txt",
						 False) # self.n_cpu)
						for i in iterator)
		else:
			for i in iterator:
				self.docker_handler.molecular_docking(
					self.docker_handler.docking_config,
					self.docker_handler.docking_receptor,
					work_path+names[i]+".pdbqt",
					docked_path+names[i]+"_docked.pdbqt",
					docked_path+names[i]+"_logs.txt",
					False) # self.n_cpu)
		print(time.time(), "Molecular docking done")
		# parse
		paths = [docked_path+names[i]+"_logs.txt" for i in range(len(names))]
		parsed_dict = pd.DataFrame(self.docker_handler.parse_vina_results(paths))
		# index
		done = []
		for i, row in parsed_dict.iterrows():
			if row["name"] not in done: 
				done.append(row["name"])
				self.index_table.loc[self.index_table["name"] == row["name"], "gen_ev"]= self.gen_counter
				self.index_table.loc[self.index_table["name"] == row["name"], "best"]  = np.min (parsed_dict[parsed_dict["name"] == row["name"]]["energy"].values)
				self.index_table.loc[self.index_table["name"] == row["name"], "mean"]  = np.mean(parsed_dict[parsed_dict["name"] == row["name"]]["energy"].values)
		# save current state to disk 
		aux_destin = docked_path+"docked_gen_"+str(self.gen_counter)+".csv"
		parsed_dict.to_csv(aux_destin, sep=",", index=False)
		aux_destin = docked_path+"table_gen_"+str(self.gen_counter)+".csv"
		self.index_table.to_csv(aux_destin, sep=",", index=False)


	def select_best(self, scale=0.3, g=None):
		''' Selects the best 50 smiles from a generation. 
			Adjusted score to boost small and diverse molecs.
			Inputs: 
			* scale: exponential weight for adjusted coefficients.
			* n: int. optional. generation ot select best molecs from. 
		'''
		if g is None: 
			g = self.gen_counter

		indexs     = self.index_table["gen_ev"] == g
		candidates = self.index_table[indexs]
		# prepare features
		mols = [Chem.MolFromSmiles(x) for x in candidates["smiles"].values]
		fingerprints     = [Chem.RDKFingerprint(mol) for mol in mols]
		# calculate adjusted coefficients
		weight_coeff     = [(900/Descriptors.MolWt(mol))**scale for mol in mols]
		similarity_coeff = []
		for i in range(len(mols)): 
			max_sim = np.max(DataStructs.BulkTanimotoSimilarity(fingerprints[i], 
					  [fingerprints[x] for x in range(len(mols)) if x != i]))
			similarity_coeff.append((1/max_sim)**scale)
		adjusted_coeff   = candidates["best"].values *\
						   np.array(weight_coeff) * np.array(similarity_coeff)
		# add overall score
		self.index_table.loc[indexs, "w_adj"] = weight_coeff
		self.index_table.loc[indexs, "s_adj"] = similarity_coeff
		self.index_table.loc[indexs, "adj"]   = adjusted_coeff
		# select best values + add to best historic
		gen_best = self.index_table[indexs].sort_values("adj")["smiles"].values
		self.historic_best += list(gen_best[:self.retain])


	def manager(self, n_gens, models_path="models/", results_path="results/",
					          work_path="aux_docking/", docked_path="docked/"):
		''' Manages the general process. '''

		print(time.time(), "Running for {0} gens".format(n_gens))
		for i in range(1+self.n_start_gen, n_gens+1):
			if not self.skip_eval:
				# evaluates a generation
				smiles, names = self.select_candidates()
				self.evaluate_pop(smiles, names, work_path, docked_path)
				# save index_table
				self.index_table.to_csv(results_path+"index_table_gen_{0}.csv"\
										.format(self.gen_counter), 
										sep=",", index=False)
			else:
				for g in range(i):
					self.select_best(g=g)
				self.skip_eval = False
			# increase number of gens
			if i < n_gens:
				self.select_best()
				save_smiles(self.historic_best,
							results_path+"historic_best{0}.txt".format(\
							self.gen_counter))
				# select best of generation
				self.gen_counter += 1
				# retrain model
				data = self.index_table.sort_values("adj")["smiles"]
				data = data.values[self.retain * self.gen_counter]
				his  = self.model_handler.retrain(epochs=5,
				 								  data=data,
				 								  params_gen=self.params_gen)
				self.model_handler.model.save(models_path+"gen_{0}.h5"\
				 							  .format(self.gen_counter))
				# new_gen
				self.new_gen()
		return True


if __name__ == "__main__":
	print("Starting the drug refinement process.")
	# get config and necessary variables
	try:
		CONFIG       = read_json(str(sys.argv[1]))
	except:
		CONFIG       = read_json("config.json")

	predictor_params = {"max_length"    : CONFIG["smiles_max_length"],
						"sampling_temp" : CONFIG["sampling_temp"], 
						"end_char"      : CONFIG["end_char"],
						"batch_predict" : CONFIG["batch_predict"],
						"verbose"       : CONFIG["sampling_verbose"]}
	model_handler    = ModelHandler(model=CONFIG["model"],
								    predictor_params=predictor_params)
	docker_handler   = DockerHandler(vina_path=CONFIG["vina_path"],
									 openbabel_path=CONFIG["openbabel_path"],
									 obminimize_path=CONFIG["obminimize_path"],
								     docking_config=CONFIG["docking_config"], 
								     docking_receptor=CONFIG["docking_receptor"])
	# params - prepare variables for optimizer
	MODEL_HANDLER  = model_handler
	DOCKER_HANDLER = docker_handler
	DATA_SOURCE    = CONFIG["data_source"]
	BASE_SIM       = CONFIG["optimizer_params"]["base_sim"]
	FIRST_SMILES   = CONFIG["optimizer_params"]["first_smiles"]
	GEN_SMILES     = CONFIG["optimizer_params"]["gen_smiles"]
	RETAIN         = CONFIG["optimizer_params"]["retain"]
	START_SMILES   = CONFIG["optimizer_params"]["start_smiles"]
	START_GEN      = CONFIG["optimizer_params"]["start_gen"]
	N_START_GEN    = CONFIG["optimizer_params"]["n_start_gen"]
	PARALLEL       = CONFIG["optimizer_params"]["parallel"]
	N_CPU          = CONFIG["optimizer_params"]["n_cpu"]
	# vars for manager
	N_GENS         = CONFIG["optimizer_params"]["n_gens"]
	# paths
	MODELS_PATH    = CONFIG["exp_name"]+CONFIG["paths"]["models_path"]
	RESULTS_PATH   = CONFIG["exp_name"]+CONFIG["paths"]["results_path"]
	WORK_PATH      = CONFIG["exp_name"]+CONFIG["paths"]["work_path"]
	DOCKED_PATH    = CONFIG["exp_name"]+CONFIG["paths"]["docked_path"]
	# instantiate parameter
	optimizer = GensOptimizer(config=CONFIG,
							  model_handler=MODEL_HANDLER,
							  docker_handler=DOCKER_HANDLER,
							  data_source=DATA_SOURCE,
							  base_sim=BASE_SIM,
							  first_smiles=FIRST_SMILES,
							  gen_smiles=GEN_SMILES,
							  retain=RETAIN,
							  start_smiles=START_SMILES,
							  start_gen=START_GEN,
							  n_start_gen=N_START_GEN,
							  parallel=PARALLEL,
							  n_cpu=N_CPU)
	optimizer.manager(N_GENS, MODELS_PATH, RESULTS_PATH, WORK_PATH, DOCKED_PATH)

print("Finished!")