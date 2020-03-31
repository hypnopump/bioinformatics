#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json

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
from smiles_predictor import *
from smiles_generator import *


class ModelHandler():
	def __init__(self, model, predictor=False, predictor_params=False):
		''' Handles generation, retraining, and saving the models.
			Inputs: 
			* model: tf.keras model instance
			* predictor: SmilesPredictor instance. 
			* predictor_params: dict. arguments for the SmilesGenerator class.
		'''
		if isinstance(model, str):
			self.model = tf.keras.models.load_model(model)
		else: 
			self.model = model
		# handles predictor class
		if not predictor and isinstance(predictor_params, dict):
			self.predictor = SmilesPredictor(self.model, **predictor_params)
		else: 
			self.predictor = predictor


	def generate(self, n=100, valid=False):
		''' Wraps the generation. 
			Inputs: 
			* n: int. number of smiles to predict for.
			* valid: boolean. whether to check for validity of predictions.
		'''
		return self.predictor.smiles_predict(num=n, valid=valid)


	def retrain(self, epochs, data, datagen=False, params_gen=False,
				val_data=None, val_datagen=None,
				compilation="categorical_crossentropy", opt="adam"):
		''' Retrains a model with the new data. 
			Inputs: 
			* epochs: int. number of epochs to train.
			* data: iterable. smiles for training
			* datagen: optional. SmilesGenerator instance.
			* params_gen: optional. dict. params for the SmilesGenerator class.
			* val_data: optional. iterable. smiles for validation
			* val_datagen: optional. SmilesGenerator instance. 
			* compilation: Boolean. Whether to compile the model before training.
			* opt: Optimizer for model compilation. Only if compilation=True.
		'''
		if data and not datagen and isinstance(params_gen, dict):
			datagen = SmilesGenerator(smiles = data, **params_gen)
		if compilation: 
			self.model.compile(loss=compilation, optimizer=opt)
		# train model
		if val_data is None:
			his = self.model.fit_generator(datagen, 
					   steps_per_epoch=int(len(data)/params_gen["batch_size"]),
					   epochs=epochs, 
					   verbose=params_gen["verbose"])
		else:
			his = self.model.fit_generator(datagen, 
					   steps_per_epoch=int(len(data)/params_gen["batch_size"]),
					   epochs=epochs, 
					   validation_data=val_datagen,
				   	   validation_steps=int(len(val_data) / \
				   	   					    params_gen["batch_size"]),
					   verbose=params_gen["verbose"])
			
		return his