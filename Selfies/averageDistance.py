from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import RDConfig
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import QED, AllChem
from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem import MolToSmiles as mol2smi
import selfies

_get_fp = lambda x: Chem.RDKFingerprint(x)

def sanitize_smiles(smi):
    '''Return a canonical smile representation of smi
    
    Parameters:
    smi (string) : smile string to be canonicalized 
    
    Returns:
    mol (rdkit.Chem.rdchem.Mol) : RdKit mol object                          (None if invalid smile string smi)
    smi_canon (string)          : Canonicalized smile representation of smi (None if invalid smile string smi)
    conversion_successful (bool): True/False to indicate if conversion was  successful 
    '''
    # try:
    mol = smi2mol(smi, sanitize=True)
    smi_canon = mol2smi(mol, isomericSmiles=False, canonical=True)
    return smi_canon
    # except:
    #     return (None, None, False)

def get_fp(mol_or_smi):
    if type(mol_or_smi) in [Chem.rdchem.Mol, Chem.rdchem.RWMol]:
        _mol = mol_or_smi
    elif type(mol_or_smi) == str:
        _mol = Chem.MolFromSmiles(mol_or_smi)
    else:
        raise ValueError("This type is not allowed.")
    return _get_fp(_mol)

def cal_avg_dist(solutions):

    dist_sum = 0
    _n = len(solutions) #solutions is pandas dataframe?

    for i in range(_n - 1):
        for j in range(i + 1, _n):
            fps1 = get_fp(solutions[i])
            fps2 = get_fp(solutions[j])
            dist = TanimotoSimilarity(fps1, fps2)
            dist_sum += dist


    return dist_sum / (_n * (_n - 1) / 2)

f = open("final_bank_0.99730_0.900.csv") #--------put csv here---------------#
Molecules = []

for row in f:
    row = row.split(',')
    Molecules.append(row[0])

Molecules.pop(0) #gets rid of first element since it's a label for us

for i in range(0, len(Molecules)):
    #-------------this line is for Selfie input----------------#
    #Molecules[i] = selfies.decoder(Molecules[i])
    Molecules[i] = sanitize_smiles(Molecules[i])

print(cal_avg_dist(Molecules))