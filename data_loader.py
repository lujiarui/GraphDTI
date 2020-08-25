"""Data Loader for DTI datasets
Kd is used as label(numerical or threshold)

- The DrugBank dataset can be found in https://www.drugbank.ca/.
- The BindingDB dataset is hosted in https://www.bindingdb.org/bind/index.jsp.
	- url = 'https://www.bindingdb.org/bind/downloads/BindingDB_All_2020m2.tsv.zip'
- The Davis Dataset can be found in http://staff.cs.utu.fi/~aatapa/data/DrugTarget/.
- The KIBA dataset can be found in https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0209-z. 
- The Drug Target Common Dataset can be found in https://drugtargetcommons.fimm.fi/.
	- url = 'https://drugtargetcommons.fimm.fi/static/Excell_files/DTC_data.csv'
- The COVID-19 Dataset including SARS-CoV, Broad Repurposing Hub can be found in https://www.aicures.mit.edu/data; and https://pubchem.ncbi.nlm.nih.gov/bioassay/1706. 

- Use some existing files from https://github.com/yangkevin2/coronavirus_data
- Use the SMILES, protein sequence from DeepDTA github repo: https://github.com/hkmztrk/DeepDTA/tree/master/data.
"""

import numpy as np
import pandas as pd
import json
import os
try:
	import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

def unit_transformer(val, from_, to_):
	"""
    Transform the unit of label(val) between nM and pK (-log10K)
    """
    # basis as nM
	if from_ == 'nM':
		val = val
	elif from_ == 'p':
		val = 10**(-val) / 1e-9
	
	if to_ == 'p':
		val = -np.log10(val*1e-9 + 1e-10)
	elif to_ == 'nM':
		val = val
	return val



# ##################################################
# ########### Common Benchmarks ####################
# ##################################################

def load_DrugBank(path='./data'):
	"""Intrinsic Binary"""
	print('Loading Dataset from ' + path + '/DrugBank')
	drug_path = path + '/DrugBank/drugs.json'
	target_path = path + '/DrugBank/targets.json'
	with open(drug_path, 'r') as f:
		drugs = json.load(f)
	with open(target_path, 'r') as f:
		targets = json.load(f)
	
	print('Done!')
	return np.array(drugs), np.array(targets)



def load_BindingDB(path='./data', val='Kd', binary=False, convert_to_log=True, threshold=30.0):
	"""Process data using pandas dataframe
	"""
	path = path + '/BindingDB/BindingDB_All.tsv'

	print('Loading Dataset from ' + path)
	df = pd.read_csv(path, sep = '\t', error_bad_lines=False)

	df = df[df['Number of Protein Chains in Target (>1 implies a multichain complex)'] == 1.0]
	df = df[df['Ligand SMILES'].notnull()]

	if val == 'Kd':
		idx_str = 'Kd (nM)'
	elif val == 'IC50':
		idx_str = 'IC50 (nM)'
	elif val == 'Ki':
		idx_str = 'Ki (nM)'
	elif val == 'EC50':
		idx_str = 'EC50 (nM)'
	else:
		print('Please select Kd, Ki, IC50 or EC50')
		return

	df_candidate = df[df[idx_str].notnull()]
	df_candidate = df_candidate[['BindingDB Reactant_set_id', 'Ligand InChI', 'Ligand SMILES',\
					  'PubChem CID', 'UniProt (SwissProt) Primary ID of Target Chain',\
					  'BindingDB Target Chain  Sequence', idx_str]]
	df_candidate.rename(columns={'BindingDB Reactant_set_id':'ID',
							'Ligand SMILES':'SMILES',
							'Ligand InChI':'InChI',
							'PubChem CID':'PubChem_ID',
							'UniProt (SwissProt) Primary ID of Target Chain':'UniProt_ID',
							'BindingDB Target Chain  Sequence': 'Target Sequence',
							idx_str: 'Label'}, 
							inplace=True)

	df_candidate['Label'] = df_candidate['Label'].str.replace('>', '')
	df_candidate['Label'] = df_candidate['Label'].str.replace('<', '')
	df_candidate['Label'] = df_candidate['Label'].astype(float)
	
	# have at least uniprot or pubchem ID
	df_candidate = df_candidate[df_candidate.PubChem_ID.notnull() | df_candidate.UniProt_ID.notnull()]
	df_candidate = df_candidate[df_candidate.InChI.notnull()]

	df_candidate = df_candidate[df_candidate.Label <= 10000000.0]
	print('There are ' + str(len(df_candidate)) + ' drug target pairs.')

	if binary:
		print('[INFO] Default binary threshold for the binding affinity scores are 30, \
			   you can adjust it by using the "threshold" parameter')
		val = [1 if i else 0 for i in df_candidate.Label.values < threshold]
	else:
		if convert_to_log:
			print('[INFO] Default set to logspace (nM -> p) for easier regression')
			val = unit_transformer(df_candidate.Label.values, 'nM', 'p') 
		else:
			val = df_candidate.Label.values

	return df_candidate.SMILES.values, df_candidate['Target Sequence'].values, np.array(val)

def load_DAVIS(path='./data', binary=False, convert_to_log=True, threshold=30.0):
	"""Process Davis database and return binary/numerical labels
	"""
	print('Processing Davis database.')

	affinity = pd.read_csv(path + '/DAVIS/affinity.txt', header=None, sep=' ')

	with open(path + '/DAVIS/target_seq.txt') as f:
		target = json.load(f)

	with open(path + '/DAVIS/SMILES.txt') as f:
		drug = json.load(f)

	# ignore the keys(sequentially arranged)
	target = list(target.values())
	drug = list(drug.values())
	
	# lists to be extracted
	SMILES = []
	targets_sequences = []
	val = []

	for i in range(len(drug)):
		for j in range(len(target)):
			SMILES.append(drug[i])
			targets_sequences.append(target[j])
			val.append(affinity.values[i, j])

	if binary:
		print('[INFO] Default binary threshold for the binding affinity scores are 30, \
			   you can adjust it by using the "threshold" parameter')
		val = [1 if i else 0 for i in np.array(val) < threshold]
	else:
		if convert_to_log:
			print('[INFO] Default set to logspace (nM -> p) for easier regression')
			val = unit_transformer(np.array(val), 'nM', 'p') 
		else:
			val = val
	print('Done!')
	return np.array(SMILES), np.array(targets_sequences), np.array(val)


def load_KIBA(path='./data', binary=False, threshold=9):
	"""Process KIBA Database and return binary/numerical labels
	"""
	print('Processing KIBA database.')

	if not os.path.exists(path):
	    os.makedirs(path)
 

	affinity = pd.read_csv(path + '/KIBA/affinity.txt', header=None, sep = '\t')
	affinity = affinity.fillna(-1)

	with open(path + '/KIBA/target_seq.txt') as f:
		target = json.load(f)

	with open(path + '/KIBA/SMILES.txt') as f:
		drug = json.load(f)
	
	# ignore the keys(sequentially arranged)
	target = list(target.values())
	drug = list(drug.values())
	
	# lists to be extracted
	SMILES = []
	targets_sequences = []
	val = []

	for i in range(len(drug)):
		for j in range(len(target)):
			if affinity.values[i, j] != -1:
				SMILES.append(drug[i])
				targets_sequences.append(target[j])
				val.append(affinity.values[i, j])

	if binary:
		print('[INFO] Note that KIBA is not suitable for binary classification as it is a modified score. \
			   Default binary threshold for the binding affinity scores are 9, \
			   you should adjust it by using the "threshold" parameter')
		val = [1 if i else 0 for i in np.array(val) < threshold]
	else:
		val = val

	print('Done!')
	return np.array(SMILES), np.array(targets_sequences), np.array(val)



