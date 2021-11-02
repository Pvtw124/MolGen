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



def crossover(selfie1, selfie2, min_contribution = 4):
    selfie1 = get_selfie_chars(selfie1)
    selfie2 = get_selfie_chars(selfie2)
    n = len(selfie1)
    m = len(selfie2)
    x = max(min_contribution, n-m)#how many chars must selfie1 use to ensure that we can create a molecule with the length of selfie1
    from_one = np.random.randint(x, n-1)#number of molecules to take from selfie1
    from_two = n-from_one #number of molecules that must be taken from selfie2
    #notice from_one + from_two = n (the lengths must always add to the length of selfie1)
    #break 2 molecules into 4 parts
    left_even = selfie1[0:from_one]
    right_even = selfie2[m-from_two:m]
    left_odd = selfie2[0:from_two]
    right_odd = selfie1[n-from_one:n]

    selfie_mutated_even = left_even + right_even
    selfie_mutated_odd = left_odd + right_odd
    selfie_mutated_even = "".join(x for x in selfie_mutated_even)
    selfie_mutated_odd = "".join(x for x in selfie_mutated_odd)
    smiles_even = selfies.decoder(selfie_mutated_even)
    smiles_odd = selfies.decoder(selfie_mutated_odd)
    mol_even, smiles_canon_even, done = sanitize_smiles(smiles_even)
    mol_odd, smiles_canon_odd, done = sanitize_smiles(smiles_odd)
    return selfie_mutated_even, smiles_canon_even, selfie_mutated_odd, smiles_canon_odd

# def crossoverRight(selfie1, selfie2, min_contribution = 4):
#     selfie1 = get_selfie_chars(selfie1)
#     selfie2 = get_selfie_chars(selfie2)
#     n = len(selfie1)
#     m = len(selfie2)
#     x1 = max(min_contribution, n-m)
#     idx1 = np.random.randint(x, n-1)

#     left2 = selfie2[0:(n-idx)]
#     right2 = selfie1[idx:n]
#     left1 = selfie1[0:idx]
#     right1 = selfie2[m-(n-idx):m]
#     selfie_mutated = left + right
#     selfie_mutated = "".join(x for x in selfie_mutated)
#     smiles = selfies.decoder(selfie_mutated)
#     mol, smiles_canon, done = sanitize_smiles(smiles)
#     return selfie_mutated, smiles_canon
# def crossover_right(selfie1, selfie2, min_contribution = 4):
#     n = len(selfie1)
#     m = len(selfie2)
#     x = max(min_contribution, n-m)
#     idx = np.random.randint(x, n-1)
#     right = selfie1[idx:n]
#     left = selfie2[m-(n-idx):m-1]
#     selfie_mutated = left + right
#     selfie_mutated = "".join(x for x in selfie_mutated)
#     smiles = selfies.decoder(selfie_mutated)
#     mol, smiles_canon, done = sanitize_smiles(smiles)
#     return selfie_mutated, smiles_canon