import selfies
import numpy as np
import rdkit
from rdkit.Chem import Draw
from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem import MolToSmiles as mol2smi
from rdkit.Chem import Descriptors

max_molecules_len = 140
alphabet = ['[Branch1_1]', '[Branch1_2]','[Branch1_3]', '[epsilon]', '[Ring1]', '[Ring2]', '[Branch2_1]', '[Branch2_2]', '[Branch2_3]', '[F]', '[O]', '[=O]', '[N]', '[=N]', '[#N]', '[C]', '[=C]', '[#C]', '[S]', '[=S]', '[C][=C][C][=C][C][=C][Ring1][Branch1_1]']

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
    return (mol, smi_canon, True)
    # except:
    #     return (None, None, False)

def get_selfie_chars(selfie):
    chars_selfie = [] # A list of all SELFIE sybols from string selfie
    while selfie != '':
        chars_selfie.append(selfie[selfie.find('['): selfie.find(']')+1])
        selfie = selfie[selfie.find(']')+1:]
    return chars_selfie


#mutations
def add(selfie):
    chars_selfie = get_selfie_chars(selfie)

    random_index = np.random.randint(len(chars_selfie)+1)
    random_character = np.random.choice(alphabet, size=1)[0]
    selfie_mutated_chars = chars_selfie[:random_index] + [random_character] + chars_selfie[random_index:]

    selfie_mutated = "".join(x for x in selfie_mutated_chars)
    sf = "".join(x for x in chars_selfie)
    # try:
    smiles = selfies.decoder(selfie_mutated)
    mol, smiles_canon, done = sanitize_smiles(smiles)
        # if len(smiles_canon) > max_molecules_len or smiles_canon=="":
        #     done = False
        # if done:
        #     valid = True
        # else:
        #     valid = False
    # except:
    #     valid=False
    #     if fail_counter > 1 and write_fail_cases == True:
    #         f = open("selfie_failure_cases.txt", "a+")
    #         f.write('Tried to mutate SELFIE: '+str(sf)+' To Obtain: '+str(selfie_mutated) + '\n')
    #         f.close()
    return (selfie_mutated, smiles_canon)

def remove(selfie):
    chars_selfie = get_selfie_chars(selfie)
    index = np.random.randint(len(chars_selfie))
    selfie_mutated_chars = chars_selfie
    selfie_mutated_chars.pop(index)
    

    selfie_mutated = "".join(x for x in selfie_mutated_chars)
    sf = "".join(x for x in chars_selfie)
    #try:
    smiles = selfies.decoder(selfie_mutated)
    mol, smiles_canon, done = sanitize_smiles(smiles)

    # if len(smiles_canon) > max_molecules_len or smiles_canon=="":
    #     done = False
    # if done:
    #     valid = True
    # else:
    #     valid = False
    # except:
    #     valid=False
    #     if fail_counter > 1 and write_fail_cases == True:
    #         f = open("selfie_failure_cases.txt", "a+")
    #         f.write('Tried to mutate SELFIE: '+str(sf)+' To Obtain: '+str(selfie_mutated) + '\n')
    #         f.close()
    return (selfie_mutated, smiles_canon)


    

def replace(selfie):
    chars_selfie = get_selfie_chars(selfie)
    random_index = np.random.randint(len(chars_selfie))
    random_character = np.random.choice(alphabet, size=1)[0]
    if random_index == 0:
        selfie_mutated_chars = [random_character] + chars_selfie[random_index+1:]
    else:
        selfie_mutated_chars = chars_selfie[:random_index] + [random_character] + chars_selfie[random_index+1:]

    selfie_mutated = "".join(x for x in selfie_mutated_chars)
    sf = "".join(x for x in chars_selfie)
    # try:
    smiles = selfies.decoder(selfie_mutated)
    mol, smiles_canon, done = sanitize_smiles(smiles)
        # if len(smiles_canon) > max_molecules_len or smiles_canon=="":
        #     done = False
        # if done:
        #     valid = True
        # else:
        #     valid = False
    # except:
    #     valid=False
    #     if fail_counter > 1 and write_fail_cases == True:
    #         f = open("selfie_failure_cases.txt", "a+")
    #         f.write('Tried to mutate SELFIE: '+str(sf)+' To Obtain: '+str(selfie_mutated) + '\n')
    #         f.close()
    return (selfie_mutated, smiles_canon)



#def crossover(selfie, selfie2):