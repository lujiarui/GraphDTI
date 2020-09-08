from rdkit import Chem
from rdkit.Chem import AllChem

def getFingerPrint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    morgan_fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits).ToBitString()
    return [smiles, morgan_fingerprint]
