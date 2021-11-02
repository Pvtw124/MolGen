from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import RDConfig
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import QED, AllChem
from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem import MolToSmiles as mol2smi

_get_fp = lambda x: Chem.RDKFingerprint(x)

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
    min_dist = 10
    max_dist = 0
    _n = len(solutions)

    for i in range(_n - 1):
        for j in range(i + 1, _n):
            fps1 = get_fp(solutions[i, 1])
            fps2 = get_fp(solutions[j, 1])
            dist = TanimotoSimilarity(fps1, fps2)
            dist_sum += dist
            if dist < min_dist:
                min_dist = dist
            if dist > max_dist:
                max_dist = dist

    return dist_sum / (_n * (_n - 1) / 2)  # , min_dist, max_dist

f = open("final_bank_0.99730_0.900.csv")

Molecules = []
for row in f:
    row = row.split(',')
    Molecules.append(row[0])

cal_avg_dist(Molecules)
