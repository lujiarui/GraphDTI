import numpy as np
import pandas as pd
import json
import os 

# ##################################################
# ########### General local loader #################
# ##################################################


def read_file_training_dataset_bioassay(path):
	"""First line has format as:
		Target_seq
	rest of lines have format as:
		SMILES score   
	""" 
	try:
		file = open(path, "r")
	except:
		print('[Error] Path Not Found!')
	target = file.readline()
	if target[-1:] == '\n':
		target = target[:-1]        
	X_drug = []
	y = []
	for aline in file:
		values = aline.split()
		X_drug.append(values[0])
		y.append(float(values[1]))
	file.close()
	return np.array(X_drug), target, np.array(y)

def read_file_training_dataset_drug_target_pairs(path):
	"""a line has format as:
		SMILES Target_seq score
	"""
	try:
		file = open(path, "r")
	except:
		print('[Error] Path Not Found!')
	X_drug = []
	X_target = []
	y = []
	for aline in file:
		values = aline.split()
		X_drug.append(values[0])
		X_target.append(values[1])
		y.append(float(values[2]))
	file.close()
	return np.array(X_drug), np.array(X_target), np.array(y)

def read_file_virtual_screening_drug_target_pairs(path):
	"""a line has format as:
		SMILES Target_seq    
	""" 
	try:
		file = open(path, "r")
	except:
		print('[Error] Path Not Found!')
	X_drug = []
	X_target = []
	for aline in file:
		values = aline.split()
		X_drug.append(values[0])
		X_target.append(values[1])
	file.close()
	return np.array(X_drug), np.array(X_target)


def read_file_repurposing_library(path):
	"""a line has format as:
		drug_names SMILES 
	""" 
	try:
		file = open(path, "r")
	except:
		print('[Error] Path Not Found!')
	X_drug = []
	X_drug_names = []
	for aline in file:
		values = aline.split()
		X_drug.append(values[1])
		X_drug_names.append(values[0])
	file.close()
	return np.array(X_drug), np.array(X_drug_names)


def read_file_target_sequence(path):
	"""a line has format as:
		target_name Target_seq    
	""" 
	try:
		file = open(path, "r")
	except:
		print('[Error] Path Not Found!')
	values = file.readline().split()
	file.close()
	return values[1], values[0]