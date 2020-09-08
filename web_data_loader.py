def load_AID1706_txt_file(path = './data'):
	print('Beginning Processing...')

	if not os.path.exists(path):
		os.makedirs(path)  
        
	url = 'https://deeppurpose.s3.amazonaws.com/AID1706.txt'
	saved_path_data = wget.download(url, path)
	return saved_path_data

def load_AID1706_SARS_CoV_3CL(path = './data', binary = True, threshold = 15, balanced = True, oversample_num = 30, seed = 1):
	print('Beginning Processing...')

	if not os.path.exists(path):
		os.makedirs(path)

	target = 'SGFKKLVSPSSAVEKCIVSVSYRGNNLNGLWLGDSIYCPRHVLGKFSGDQWGDVLNLANNHEFEVVTQNGVTLNVVSRRLKGAVLILQTAVANAETPKYKFVKANCGDSFTIACSYGGTVIGLYPVTMRSNGTIRASFLAGACGSVGFNIEKGVVNFFYMHHLELPNALHTGTDLMGEFYGGYVDEEVAQRVPPDNLVTNNIVAWLYAAIISVKESSFSQPKWLESTTVSIEDYNRWASDNGFTPFSTSTAITKLSAITGVDVCKLLRTIMVKSAQWGSDPILGQYNFEDELTPESVFNQVGGVRLQ'
	url = 'https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?query=download&record_type=datatable&actvty=all&response_type=save&aid=1706'
	saved_path_data = wget.download(url, path)

	url = 'https://drive.google.com/uc?export=download&id=1eipPaFrg-mVULoBhyp2kvEemi2WhDxsM'
	saved_path_conversion = wget.download(url, path)

	df_data = pd.read_csv(saved_path_data)
	df_conversion = pd.read_csv(saved_path_conversion)
	val = df_data.iloc[4:][['PUBCHEM_CID','PUBCHEM_ACTIVITY_SCORE']]

	val['binary_label'] = 0
	val['binary_label'][(val.PUBCHEM_ACTIVITY_SCORE >= threshold) & (val.PUBCHEM_ACTIVITY_SCORE <=100)] = 1

	if balanced:
		val = pd.concat([val[val.binary_label==0].sample(n = len(val[val.binary_label==1]) * oversample_num, replace = False, random_state = seed), pd.concat([val[val.binary_label==1]]*oversample_num, ignore_index=True)]).sample(frac = 1, replace = False, random_state = seed).reset_index(drop = True)

	cid2smiles = dict(zip(df_conversion[['cid','smiles']].values[:, 0], df_conversion[['cid','smiles']].values[:, 1]))
	X_drug = [cid2smiles[i] for i in val.PUBCHEM_CID.values]    

	if binary:
		print('Default binary threshold for the binding affinity scores is 15, recommended by the investigator')
		y = val.binary_label.values
	else:
		y = val.PUBCHEM_ACTIVITY_SCORE.values

	print('Done!')
	return np.array(X_drug), target, np.array(y)

def load_broad_repurposing_hub(path = './data'):
	url = 'https://deeppurpose.s3.amazonaws.com/broad.csv'
	if not os.path.exists(path):
	    os.makedirs(path)
	saved_path_data = wget.download(url, path)
	df = pd.read_csv(saved_path_data)
	df = df.fillna('UNK')
	return df.smiles.values, df.title.values, df.cid.values.astype(str)

def load_antiviral_drugs(path = './data', no_cid = False):
	url = 'https://deeppurpose.s3.amazonaws.com/antiviral_drugs.csv'
	if not os.path.exists(path):
		os.mkdir(path)
	saved_path_data = wget.download(url, path)
	df = pd.read_csv(saved_path_data)
	if no_cid:
		return df.SMILES.values, df[' Name'].values
	else:
		return df.SMILES.values, df[' Name'].values, df['Pubchem CID'].values

def load_IC50_Not_Pretrained(path = './data', n=500): 
	print('Downloading...')    
	url = 'https://deeppurpose.s3.amazonaws.com/IC50_not_Kd.csv'
	if not os.path.exists(path):
	    os.makedirs(path)
	saved_path_data = wget.download(url, path)
	df = pd.read_csv(saved_path_data).sample(n = n, replace = False).reset_index(drop = True)
	return df['Target Sequence'].values, df['SMILES'].values

def load_IC50_1000_Samples(path = './data', n=100): 
	print('Downloading...')    
	url = 'https://deeppurpose.s3.amazonaws.com/IC50_samples.csv'
	if not os.path.exists(path):
	    os.makedirs(path)
	saved_path_data = wget.download(url, path)
	df = pd.read_csv(saved_path_data).sample(n = n, replace = False).reset_index(drop = True)
	return df['Target Sequence'].values, df['SMILES'].values

