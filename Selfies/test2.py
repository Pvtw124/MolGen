import selfies
import numpy as np



def get_selfie_chars(selfie):
    '''Obtain a list of all selfie characters in string selfie
    
    Parameters: 
    selfie (string) : A selfie string - representing a molecule 
    
    Example: 
    get_selfie_chars('[C][=C][C][=C][C][=C][Ring1][Branch1_1]')
    ['[C]', '[=C]', '[C]', '[=C]', '[C]', '[=C]', '[Ring1]', '[Branch1_1]']
    
    Returns:
    chars_selfie: list of selfie characters present in molecule selfie
    '''
    chars_selfie = [] # A list of all SELFIE sybols from string selfie
    while selfie != '':
        chars_selfie.append(selfie[selfie.find('['): selfie.find(']')+1])
        selfie = selfie[selfie.find(']')+1:]
    return chars_selfie


smiles = "Cc1cc(C(=O)N[C@H]2CCCN(Cc3ocnc3C)[C@H]2C)on1"
selfie = selfies.encoder(smiles)

selfie_chars = get_selfie_chars(selfie)
n = len(selfie_chars)
index = np.random.randint(n)
print("popped")
popped = selfie_chars.pop(index)
print(popped)
modifiedSelfie = "".join(x for x in selfie_chars)


modifiedSmi = selfies.decoder(modifiedSelfie)

print(f'original Smiles {smiles}')
print(f'original Selfies {selfie}')
print(f'New Selfies {modifiedSelfie}')
print(f'New Smiles {modifiedSmi}')