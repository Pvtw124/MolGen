import selfies
import modSELFIES
import numpy as np



# def crossover(selfie1, selfie2, min_contribution = 4):
#     n = len(selfie1)
#     m = len(selfie2)
#     x = max(min_contribution, n-m)
#     idx = np.random.randint(x, n-1)
#     left = selfie1[0:idx]
#     right = selfie2[m-(n-idx):m-1]
#     crossover = left + right
#     left = "".join(x for x in left)
#     right = "".join(x for x in right)
#     print(f"m = {m}, n = {n}, idx = {idx}")
#     print(left)
#     print(right)
#     crossover = "".join(x for x in crossover)
#     print("crossover: ", crossover)

#     return crossover

def crossover(selfie1, selfie2, min_contribution = 4):
    n = len(selfie1)
    m = len(selfie2)
    x = max(min_contribution, n-m)
    idx = np.random.randint(x, n-1)
    left = selfie1[0:idx]
    right = selfie2[m-(n-idx):m-1]
    selfie_mutated = left + right
    selfie_mutated = "".join(x for x in selfie_mutated)
    smiles = selfies.decoder(selfie_mutated)
    print(smiles)
    mol, smiles_canon, done = modSELFIES.sanitize_smiles(smiles)
    print(smiles_canon)
    return selfie_mutated, smiles_canon

def get_selfie_chars(selfie):
    chars_selfie = [] # A list of all SELFIE sybols from string selfie
    while selfie != '':
        chars_selfie.append(selfie[selfie.find('['): selfie.find(']')+1])
        selfie = selfie[selfie.find(']')+1:]
    return chars_selfie

selfie1 = '[C][C][C][=C][C][=C][C][=C][Ring1][Branch1_2][C][C][Branch1_2][C][=O][N][Branch1_1][Ring1][C][C][C][C][N][C][Branch1_2][C][=O][C][=C][C][Branch1_1][C][C][=C][O][Ring1][Branch1_2]'
selfie2 = '[C][C][=N][N][Branch1_1][Branch2_3][C][C][Branch1_1][C][F][Branch1_1][C][F][F][C][=C][Ring1][Branch2_3][C][Branch1_2][C][=O][O][C][C@@Hexpl][C][C][C][N][C][Ring1][Branch1_2][=O]'

selfie1 = get_selfie_chars(selfie1)
selfie2 = get_selfie_chars(selfie2)

crossover(selfie1, selfie2)

