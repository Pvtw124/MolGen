#!/usr/bin/env python3

import os
import sys
import time
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import selfies 

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import RDConfig
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import QED, AllChem
from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem import MolToSmiles as mol2smi
#from rdkit.Chem.rdchem import KekulizeException

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score')) #adds SA_Score to some file?

import sascorer #synthetic accessibility scorer

import ModSMI #other python file

parser = argparse.ArgumentParser(
    description="This python script is made by Yongbeom Kwon")
parser.add_argument("-i", "--input", type=str, required=True, help="")
parser.add_argument("-r",
                    "--random-seed",
                    metavar="N",
                    type=int,
                    default=None,
                    help="Random seed")
parser.add_argument(
    "--bank-size",
    metavar="N",
    type=int,
    default=100,
    help="",
)
parser.add_argument(
    "--seed-size",
    metavar="N",
    type=int,
    default=60,
    help="",
)
parser.add_argument(
    "--max-round",
    metavar="N",
    type=int,
    default=150,
    help="",
)
parser.add_argument(
    "-cvg",
    "--convergent-round",
    metavar="N",
    default=150,
    type=int,
    help=
    "Convergent round; It determines when D_cut reaches a minimum value. And also It decides diversity of molecules",
)
parser.add_argument(
    "-c",
    "--coefficient",
    metavar="Float",
    type=float,
    default=0.9,
    help="coefficient of reward function.",
)
parser.add_argument(
    "-dist",
    "--dist-coef",
    metavar="coef. of distance",
    type=float,
    default=0.90,
    help="Control Dcut",
)
parser.add_argument(
    "--target",
    metavar="SMILES",
    type=str,
    default=None,
    help="target_moleclue SMILES",
)
parser.add_argument(
    "-fp",
    "--fp-method",
    type=str,
    default="rdkit",
    help="Select Fingerprint Method (rdkit/morgan)",
)
# parser.add_argument(
#     "-nf", "--nfeatures", metavar="N", type=int, default=2, help="a number of features"
# )
parser.add_argument("-v",
                    "--verbosity",
                    action="count",
                    default=0,
                    help="print error")
args = parser.parse_args()

if args.verbosity == 0:
    rdBase.DisableLog('rdApp.*')

fp_method = args.fp_method #choose what finger print algorithm you want (for calculating where it is in chemical space)

if fp_method == "rdkit":
    _get_fp = lambda x: Chem.RDKFingerprint(x)
elif fp_method == "morgan":
    _get_fp = lambda x: AllChem.GetMorganFingerprint(x, 2)
elif fp_method == "morgan1024":
    _get_fp = lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024)
elif fp_method == "morgan2048":
    _get_fp = lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2)


def get_fp(mol_or_smi):
    if type(mol_or_smi) in [Chem.rdchem.Mol, Chem.rdchem.RWMol]:
        _mol = mol_or_smi
    elif type(mol_or_smi) == str:
        _mol = selfies.decoder(mol_or_smi)#convert to smiles
        _mol = Chem.MolFromSmiles(mol_or_smi) #change this so it converts from SElFIES
    else:
        raise ValueError("This type is not allowed.")
    return _get_fp(_mol)


# -------------------------------------------------------- #
"""
@@FEATURES
Features must return a number.
1. Set nfeatures
2. Add your functions
"""

# The number of features
nfeatures = 2

# Set your functions.
feature1 = lambda x: sascorer.calculateScore(x)  # SAS (synthetic accessibility)
feature1a = lambda x: 1 - (sascorer.calculateScore(x) - 1) / 9  # SAS (synthetic accessibility)
feature2 = lambda x: QED.default(x)  # QED (Quantative estimate of drug likeness) (QED instead of log?)
feature3 = lambda x: TanimotoSimilarity(_get_fp(x), target_fps)  # similarity https://docs.eyesopen.com/toolkits/python/graphsimtk/measure.html


def feature5(x):  # RingPenalty
    ring_list = x.GetRingInfo().AtomRings()
    if ring_list:
        ring_length = max([len(j) for j in ring_list])
        if ring_length <= 6:
            return 0
        else:
            return ring_length - 6
    else:
        return 0


# SET REWARD function #################################### #
"""
@@REWARD
Modify the cal_features function; Index 0, 1, and 2 are fixed.
You need to put the desired features after index 3 of the returned list.
Then modify obj_fn and column_name to suit your purpose.
"""
if args.target:  # TARGET Molecule VESRSION

    # Modify the cal_features function; Index 0, 1, and 2 are fixed.
    # You need to put the desired features after index 3 of the returned list.
    def cal_features(_smi, _mol):  # set data_column
        return [_smi, _mol, True, feature2(_mol), feature3(_mol)]

    # Then modify obj_fn and column_name to suit your purpose.
    target_fps = get_fp(args.target)
    sim_coef = args.coefficient
    qed_coef = 1 - sim_coef
    obj_fn = lambda x: sim_coef * x[:, 4] + qed_coef * x[:, 3]

    #            SMILES / feature1 / feature2 / obj_fn
    column_name = ["SMILES", "QED", "SIM", "TARGET"]

else:  # Non-TARGET Molecule VESRSION

    # Modify the cal_features function; Index 0, 1, and 2 are fixed.
    # You need to put the desired features after index 3 of the returned list.
    def cal_features(_smi, _mol):  # set data_column
        return [_smi, _mol, True, feature1(_mol), feature2(_mol)]

    # Then modify obj_fn and column_name to suit your purpose.
    qed_coef = args.coefficient
    sas_coef = 1 - qed_coef
    obj_fn = lambda x: qed_coef * x[:, 4] - sas_coef * x[:, 3]

    #            SMILES / feature1 / feature2 / obj_fn
    column_name = ["SMILES", "SAS", "QED", "TARGET"]

# -------------------------------------------------------- #


class ChkTime:
    def __init__(self):
        self.t = time.time()

    def get(self):
        return (time.time() - self.t) / 60


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


def cal_rnd_avg_dist(solutions, nrnd=400000): #used if bank is over 600, otherwise not random

    dist_sum = 0
    min_dist = 10
    max_dist = 0
    tmp_chk = 0

    start_chk = time.time()
    for _ in range(nrnd):

        if _ == 0:
            tmp_chk = start_chk

        mol1, mol2 = np.random.choice(solutions[:, 1], size=2)
        fps1 = get_fp(mol1)
        fps2 = get_fp(mol2)
        dist = TanimotoSimilarity(fps1, fps2)
        dist_sum += dist

        if dist < min_dist:
            min_dist = dist
        if dist > max_dist:
            max_dist = dist
        if _ % int(nrnd / 10) == 0:
            print(
                f"{_/nrnd*100:.1f}% complete {(time.time() - tmp_chk)/60} min/10%\r"
            )
            tmp_chk = time.time()

    print(f"calc. Dist total {(time.time() - start_chk)/60} min")

    return dist_sum / nrnd  # , min_dist, max_dist


def cal_array_dist(solutions1, solutions2): #not used
    """
    numpy
    :param solutions1:
    :param solutions2:
    :return:
    """

    n1 = len(solutions1)
    n2 = len(solutions2)
    dist = np.zeros([n1, n2])

    for i in range(n1):
        for j in range(n2):
            fps1 = get_fp(solutions[i, 1])
            fps2 = get_fp(solutions[j, 1])
            dist[n1, n2] = TanimotoSimilarity(fps1, fps2)

    return dist


def init_bank(file_name, nbank=None, nsmiles=None, rseed=None):

    np.random.seed(rseed)

    _df = pd.read_csv(file_name)
    df = _df[:nsmiles].values

    shuffled_index = np.random.permutation(len(df))
    tmp_bank = df[shuffled_index][:nbank]
    df = pd.DataFrame(tmp_bank, columns=_df.columns.tolist())
    df.to_csv(f"init_bank_{R_d:.5f}_{args.coefficient:.3f}.csv", index=False)

    del df, _df

    bank = np.empty([tmp_bank.shape[0], 3 + nfeatures], dtype=object)
    bank[:, 0] = tmp_bank[:, 0]  # SMILES
    bank[:, 2] = True  # usable label True

    for i, j in enumerate(bank[:, 0]):
        # mol = selfies.decoder(j)
        
        mol = Chem.MolFromSmiles(j)   #---------------------------------------------
        bank[i, 1] = mol
        Chem.Kekulize(mol)
        bank[i, 0] = Chem.MolToSmiles(mol,
                                      kekuleSmiles=True,
                                      isomericSmiles=False)
        bank[i, 3:] = cal_features(j, mol)[3:]
    # print(bank)
    my_df = pd.DataFrame(bank)
    # my_df.to_excel("bank.xlsx")
    return bank


def prepare_seed(solutions, seed):

    solutions = solutions[np.where(solutions[:, 2] == True)]  # 사용 가능한 것이 True
    x = np.argsort(obj_fn(solutions))
    solutions = x[::-1]

    # shuffled_index = np.random.permutation(len(true_solutions))
    # true_solutions = true_solutions[shuffled_index]

    if len(solutions) > nseed:
        i = 0
        if len(seed) is 0:  # First selection,
            bank[solutions[0], 2] = False
            seed.append(bank[solutions[0]])
            i += 1

        if len(solutions) < len(seed):
            print(
                f"## Solutions is less than seeds / round {round_} / iter {niter}"
            )
            raise ValueError

        while len(seed) < nseed:
            if len(solutions) == i + 1:
                print(
                    f"## Solutions is empty state / unused > nseed / # of seed: {len(seed)} / round {round_} / iter {niter}"
                )
                break
                # raise ValueError

            # dist = cal_array_dist(solutions[i, 1], seed)
            dist = np.zeros([len(seed)])
            for j in range(len(seed)):
                fps1 = get_fp(bank[solutions[i], 1])
                fps2 = get_fp(seed[j][1])
                dist[j] = TanimotoSimilarity(fps1, fps2)
            if np.max(dist) > (1 - davg):
                i += 1
                # print(f'Dcut !!! {np.max(dist)}')
                continue
            else:
                bank[solutions[i], 2] = False
                seed.append(bank[solutions[i]])
                i += 1
    else:
        for i in bank[solutions]:
            seed.append(i)
        bank[solutions[:], 2] = False
        rnd_number = np.random.permutation(len(bank))
        i = 0
        while len(seed) <= nseed:
            if len(rnd_number) == i + 1:
                print(
                    f"## Solutions is empty state / unused < nseed / # of seed: {len(seed)} / round {round_} / iter {niter}"
                )
                break
            dist = np.zeros([len(seed)])
            for j in range(len(seed)):
                fps1 = get_fp(bank[rnd_number[i], 1])
                fps2 = get_fp(seed[j][1])
                dist[j] = TanimotoSimilarity(fps1, fps2)
            if np.max(dist) > (1 - davg):
                i += 1
                # print(f'Dcut !!! {np.max(dist)}')
                continue
            else:
                seed.append(bank[rnd_number[i]])
                i += 1

    print(f"@ prepare_seed finished!")
    return np.asarray(seed)


def cut_smi(smi1, smi2, func, ring_bool):

    l_smi = None
    r_smi = None

    try:
        l_smi = func(smi1, "L", ring_bool, 4)
        r_smi = func(smi2, "R", ring_bool, 4)
    except (IndexError, ValueError):
        fail_f.write(f"{l_smi},{r_smi},piece\n")
        raise PermissionError

    return l_smi, r_smi


def get_sliced_smiles(smi1, smi2, func, ring_bool):
    l_smi = None
    r_smi = None

    gate = 0
    while not (l_smi and r_smi):
        gate += 1
        if gate > 10:
            # fail_f.write(f"{l_smi},{r_smi},np\n")
            raise PermissionError
        try:
            l_smi, r_smi = cut_smi(smi1, smi2, func, ring_bool)
        except:
            pass
    return l_smi, r_smi


def crossover_smiles(smi1, smi2, func, ring_bool):
    new_smi = None
    mol = None

    l_smi, r_smi = get_sliced_smiles(smi1, smi2, func, ring_bool)

    gate = 0
    while not mol:
        gate += 1
        if gate > 5:
            break
        try:
            new_smi = ModSMI.tight_rm_branch(l_smi, r_smi)
        except ValueError:
            continue
        mol = selfies.decoder(new_smi)
        mol = Chem.MolFromSmiles(mol) #------------------------------------

    if not mol:
        l_smi, r_smi = get_sliced_smiles(smi2, smi1, func, ring_bool)

        gate = 0
        while not mol:
            gate += 1
            if gate > 5:
                break
            try:
                new_smi = ModSMI.tight_rm_branch(l_smi, r_smi)
            except ValueError:
                continue
            
            mol = Chem.MolFromSmiles(new_smi) #-------------------------------------------------

    return new_smi, mol


def append_seed(_smi, _mol, update_solution):
    try:
        update_solution.append(cal_features(_smi, _mol))
        return 1
    except Chem.rdchem.MolSanitizeException:
        print(f"#### QED error {new_smi}")
        return 0

def sanitize_smiles(smi):
    '''Return a canonical smile representation of smi
    
    Parameters:
    smi (string) : smile string to be canonicalized 
    
    Returns:
    mol (rdkit.Chem.rdchem.Mol) : RdKit mol object                          (None if invalid smile string smi)
    smi_canon (string)          : Canonicalized smile representation of smi (None if invalid smile string smi)
    conversion_successful (bool): True/False to indicate if conversion was  successful 
    '''
    try:
        mol = smi2mol(smi, sanitize=True)
        smi_canon = mol2smi(mol, isomericSmiles=False, canonical=True)
        return (mol, smi_canon, True)
    except:
        return (None, None, False)

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


def mutations_random_grin(selfie, max_molecules_len, write_fail_cases=False):
    '''Return a mutated selfie string
    
    Mutations are done until a valid molecule is obtained 
    Rules of mutation: With a 50% propbabily, either: 
        1. Add a random SELFIE character in the string
        2. Replace a random SELFIE character with another
    
    Parameters:
    selfie            (string)  : SELFIE string to be mutated 
    max_molecules_len (int)     : Mutations of SELFIE string are allowed up to this length
    write_fail_cases  (bool)    : If true, failed mutations are recorded in "selfie_failure_cases.txt"
    
    Returns:
    selfie_mutated    (string)  : Mutated SELFIE string
    smiles_canon      (string)  : canonical smile of mutated SELFIE string
    '''
    valid=False
    fail_counter = 0
    chars_selfie = get_selfie_chars(selfie)
    while not valid:
        fail_counter += 1
                
        alphabet = ['[Branch1_1]', '[Branch1_2]','[Branch1_3]', '[epsilon]', '[Ring1]', '[Ring2]', '[Branch2_1]', '[Branch2_2]', '[Branch2_3]', '[F]', '[O]', '[=O]', '[N]', '[=N]', '[#N]', '[C]', '[=C]', '[#C]', '[S]', '[=S]', '[C][=C][C][=C][C][=C][Ring1][Branch1_1]']

        # Insert a character in a Random Location
        if np.random.random() < 0.5: 
            random_index = np.random.randint(len(chars_selfie)+1)
            random_character = np.random.choice(alphabet, size=1)[0]

            selfie_mutated_chars = chars_selfie[:random_index] + [random_character] + chars_selfie[random_index:]

        # Replace a random character 
        else:                      
            random_index = np.random.randint(len(chars_selfie))
            random_character = np.random.choice(alphabet, size=1)[0]
            if random_index == 0:
                selfie_mutated_chars = [random_character] + chars_selfie[random_index+1:]
            else:
                selfie_mutated_chars = chars_selfie[:random_index] + [random_character] + chars_selfie[random_index+1:]
    
        selfie_mutated = "".join(x for x in selfie_mutated_chars)
        sf = "".join(x for x in chars_selfie)
        try:
            smiles = selfies.decoder(selfie_mutated)
            mol, smiles_canon, done = sanitize_smiles(smiles)
            if len(smiles_canon) > max_molecules_len or smiles_canon=="":
                done = False
            if done:
                valid = True
            else:
                valid = False
        except:
            valid=False
            if fail_counter > 1 and write_fail_cases == True:
                f = open("selfie_failure_cases.txt", "a+")
                f.write('Tried to mutate SELFIE: '+str(sf)+' To Obtain: '+str(selfie_mutated) + '\n')
                f.close()
    return (selfie_mutated, smiles_canon)



#-----------------------------add mutations_random_grin--------------------------
def prepare_child(seed):
    #convert smiles to selfies

    # for i in range(seed.shape[0]):
    #     try:
    #         smi1 = seed[i, 0]
    update_solution = []
    for i in range(seed.shape[0]):
        for j in range(0, 40):
            selfie_input = selfies.encoder(seed[i, 0])
            _, mutated_smi = mutations_random_grin(selfie_input, 140)
            mutated_mol = Chem.MolFromSmiles(mutated_smi)
            append_seed(mutated_smi, mutated_mol, update_solution)
        
    return np.asarray(update_solution)

    # update_solution = array of cal_features(_smi, _mol)
    # 		-cal_features is just [_smi, _mol, True, QED, tanimotoSimilarity]

    #mutations_grin gives us smiles_canon and selfies
    #add to array, update_solution




    # update_solution = []
    # # print(seed[:, 0])
    # for i in range(seed.shape[0]):
    #     try:
    #         smi1 = seed[i, 0]
    #     except IndexError:
    #         print(f"seed_shape: {seed.shape} / #: {i}")
    #         raise Exception

    #     # CROSSOVER1 ###
    #     q = 0
    #     j = 0
    #     while j < nCross1:
    #         if q == nCross1 * 10:
    #             print(f"#### have problems in updating solutions @{smi1}")
    #             break
    #         try:
    #             w = np.random.randint(len(bank))
    #             smi2 = bank[w, 0]

    #             if np.random.random() >= 0.5:
    #                 new_smi, mol = crossover_smiles(
    #                     smi1, smi2, ModSMI.prepare_rigid_crossover, True)
    #             else:
    #                 new_smi, mol = crossover_smiles(
    #                     smi2, smi1, ModSMI.prepare_rigid_crossover, True)

    #             if mol:
    #                 j += append_seed(new_smi, mol, update_solution)
    #                 bank[w, 2] = False

    #         except PermissionError:
    #             q += 1

    #     # CROSSOVER2 ###
    #     q = 0
    #     j = 0
    #     while j < nCross2:
    #         if q == nCross2 * 10:
    #             print(f"#### have problems in updating solutions @{smi1}")
    #             break
    #         try:
    #             w = np.random.randint(len(bank))
    #             smi2 = bank[w, 0]

    #             if np.random.random() >= 0.5:
    #                 new_smi, mol = crossover_smiles(
    #                     smi1, smi2, ModSMI.prepare_rigid_crossover, False)
    #             else:
    #                 new_smi, mol = crossover_smiles(
    #                     smi2, smi1, ModSMI.prepare_rigid_crossover, False)

    #             if mol:
    #                 j += append_seed(new_smi, mol, update_solution)
    #                 bank[w, 2] = False

    #         except PermissionError:
    #             q += 1

    #     # REPLACE ###
    #     q = 0
    #     j = 0
    #     while j < nReplace:
    #         if q == nReplace * 10:
    #             print(f"#### have problems in updating solutions @{smi1}")
    #             break
    #         try:
    #             new_smi, mol = ModSMI.replace_atom(smi1)

    #             if mol:
    #                 j += append_seed(new_smi, mol, update_solution)

    #         except (PermissionError, Chem.rdchem.KekulizeException):
    #             q += 1

    #     # print(f"### seed_len2: {len(update_solution)}")

    #     # ADD ###
    #     q = 0
    #     j = 0
    #     while j < nAdd:
    #         if q == nAdd * 10:
    #             print(f"#### have problems in updating solutions @{smi1}")
    #             break
    #         try:
    #             new_smi, mol = ModSMI.add_atom(smi1)

    #             if mol:
    #                 j += append_seed(new_smi, mol, update_solution)

    #         except PermissionError:
    #             q += 1

    #     # print(f"### seed_len3: {len(update_solution)}")

    #     # REMOVE ###
    #     q = 0
    #     j = 0
    #     while j < nRemove:
    #         if q == nRemove * 10:
    #             print(f"#### have problems in updating solutions @{smi1}")
    #             break
    #         try:
    #             new_smi, mol = ModSMI.delete_atom(smi1)

    #             if mol:
    #                 j += append_seed(new_smi, mol, update_solution)
    #         except PermissionError:
    #             q += 1

    #return np.asarray(update_solution)


def prepare_local_child(_smi, nReplace=10, nAdd=10, nRemove=10):
    _mol = selfies.decoder(_smi)
    _mol = Chem.MolFromSmiles(mol)#---------------------------------------------------
    update_solution = [cal_features(_smi, _mol)]

    # REPLACE ###
    q = 0
    j = 0
    while j < nReplace:
        if q == nReplace * 10:
            # print(f"#### have problems in updating solutions @{_smi}")
            break
        try:
            new_smi, mol = ModSMI.replace_atom(_smi)
            if mol:
                j += append_seed(new_smi, mol, update_solution)
        except (PermissionError, Chem.rdchem.KekulizeException):
            q += 1

    # print(f"### seed_len2: {len(update_solution)}")

    # ADD ###
    q = 0
    j = 0
    while j < nAdd:
        if q == nAdd * 10:
            print(f"#### have problems in updating solutions @{_smi}")
            break
        try:
            new_smi, mol = ModSMI.add_atom(_smi)
            if mol:
                j += append_seed(new_smi, mol, update_solution)
        except PermissionError:
            q += 1

    # print(f"### seed_len3: {len(update_solution)}")

    # REMOVE ###
    q = 0
    j = 0
    while j < nRemove:
        if q == nRemove * 10:
            print(f"#### have problems in updating solutions @{_smi}")
            break
        try:
            new_smi, mol = ModSMI.delete_atom(_smi)
            if mol:
                j += append_seed(new_smi, mol, update_solution)
        except PermissionError:
            q += 1

    return np.asarray(update_solution)


def update_bank(child_solutions, local_opt=False):
    cnt_replace = 0
    bank_min = np.min(obj_fn(bank))
    child_solutions = child_solutions[obj_fn(child_solutions) > bank_min]

    if len(child_solutions) == 0:
        raise PermissionError("child solutions 가 없습니다 !")

    for i in range(len(child_solutions)):
        if local_opt:
            local_solutions = prepare_local_child(child_solutions[i, 0])
            x = np.argmax(obj_fn(local_solutions))
            fps1 = get_fp(local_solutions[x, 1])
        else:
            fps1 = get_fp(child_solutions[i, 1])

        max_similarity = 0
        max_n = None
        for _ in range(len(bank)):
            fps2 = get_fp(bank[_, 1])
            dist = TanimotoSimilarity(fps1, fps2)
            if dist > max_similarity:
                max_similarity = dist
                max_n = _

        if local_opt:
            if (1 - max_similarity) < dcut:
                if obj_fn(local_solutions[x:x + 1]) > obj_fn(
                        bank[max_n:max_n + 1]):
                    bank[max_n] = local_solutions[x:x + 1]
                    cnt_replace += 1
            else:
                _min = np.argmin(obj_fn(bank))
                if (max_similarity < 0.98) and (obj_fn(bank[_min:_min + 1]) <
                                                final_avg.mean()):
                    if obj_fn(local_solutions[x:x + 1]) > obj_fn(
                            bank[_min:_min + 1]):
                        bank[_min] = local_solutions[x:x + 1]
                        cnt_replace += 1
        else:
            if (1 - max_similarity) < dcut:
                if obj_fn(child_solutions[i:i + 1]) > obj_fn(
                        bank[max_n:max_n + 1]):
                    bank[max_n] = child_solutions[i]
                    cnt_replace += 1
            else:
                _min = np.argmin(obj_fn(bank))
                if (max_similarity < 0.98) and (obj_fn(bank[_min:_min + 1]) <
                                                final_avg.mean()):
                    if obj_fn(child_solutions[i:i + 1]) > obj_fn(
                            bank[_min:_min + 1]):
                        bank[_min] = child_solutions[i]
                        cnt_replace += 1

    return cnt_replace, len(child_solutions)


if __name__ == "__main__":

    target_value = 3
    target_round = args.convergent_round

    R_d = 10**(np.log10(2 / target_value) / int(target_round))

    nbank = args.bank_size  # number of bank conformations
    nseed = args.seed_size  # number of seeds(mating) per iteration
    max_repeat = args.max_round

    total_time = ChkTime()
    chk_load = ChkTime()

    bank = init_bank(
        args.input,
        nbank,
        rseed=args.random_seed,
    )  # zinc

    chk_load = chk_load.get()

    first_bank = bank

    origin_avg = obj_fn(bank)

    plot_list = []

    fail_f = open(f"fail_smiles.txt", "w")

    chk_calc = ChkTime()

    if nbank > 600:
        davg = cal_rnd_avg_dist(first_bank)
    else:
        davg = cal_avg_dist(first_bank)

    davg = 1 - davg
    davg = davg * args.dist_coef
    dcut = davg / 2

    final_avg = origin_avg

    chk_calc = chk_calc.get()

    with open(f"iteration.log", "w") as log_f2:
        log_f2.write(f"load_time: {chk_load:.3f} min\n")
        log_f2.write(f"dist_time: {chk_calc:.1f} min\n")
        log_f2.write(
            f"round  iter  unused  time_seed  time_child  time_update  n_replace  n_eval\n"
        )

    with open(f"message.log", "w") as log_f:
        log_f.write(f"nbank: {nbank}\n")
        log_f.write(f"nseed: {nseed}\n")
        log_f.write(f"max_repeat: {max_repeat}\n")
        log_f.write(f"R_d: {R_d:.6f} (convergent_round: {target_round})\n")
        log_f.write(f"D_avg: {davg:.3f} (similarity: {1-davg:.3f})\n")
        tmp_str = ""
        for i, j in enumerate(column_name[1:-1]):
            tmp_str += f"{j}: {first_bank[:, i+3].mean():.3f}, "
        log_f.write(f"init_bank_avg - {tmp_str[:-2]}\n")
        log_f.write(
            f"round   dcut  n_iter  obj_avg  obj_min  obj_max  n_replace  min/round\n"
        )

    save_bank = np.empty([max_repeat, bank.shape[0], bank.shape[1] - 1],
                         dtype=object)

    for round_ in range(max_repeat):
        if (round_ != 0) and (dcut > davg / 3):
            dcut *= R_d

        timechk = time.time()

        # print(f'dcut: {dcut}, davg: {davg}')
        niter = 0
        n_replace = 0
        iter_gate = True
        while iter_gate:
            seed = []
            # log_f.write(f'## SEED #### @ {np.count_nonzero(bank[:, 2] == True)} #############\n')
            time_seed = ChkTime()
            seed = prepare_seed(bank, seed)
            time_seed = time_seed.get()
            # log_f.write(f'## SEED #### @ {np.count_nonzero(bank[:, 2] == True)} #### AFTER ##\n')

            # log_f.write(f'## CHILD ### @ {np.count_nonzero(bank[:, 2] == True)} #############\n')
            time_child = ChkTime()
            child_solutions = prepare_child(seed)
            shuffled_index_ = np.random.permutation(
                child_solutions.shape[0])  # @4 에서 추가 됨.
            child_solutions = child_solutions[shuffled_index_]
            time_child = time_child.get()
            # log_f.write(f'## CHILD ### @ {np.count_nonzero(bank[:, 2] == True)} #### AFTER ##\n')

            time_update = ChkTime()
            try:
                # log_f.write(f'## BANK #### @ {np.count_nonzero(bank[:, 2] == True)} #############\n')
                # n_replace += update_bank(child_solutions, True)  # local update
                _n_replace, n_eval = update_bank(
                    child_solutions)  # non-local update
                n_replace += _n_replace
                # log_f.write(f'## BANK #### @ {np.count_nonzero(bank[:, 2] == True)} #### AFTER ##\n')
            except PermissionError:
                break
            time_update = time_update.get()
            niter += 1

            if np.count_nonzero(bank[:, 2] == True) < (nbank - nseed * 0.9):
                iter_gate = False
            with open(f"iteration.log", "a") as log_f2:
                log_f2.write(
                    f"{round_:>4}  {niter:4}  {np.count_nonzero(bank[:, 2] == True):>6}     {time_seed:6.1f}"
                    f"      {time_child:6.1f}       {time_update:6.1f}    {n_replace:7}     {n_eval}\n"
                )

        final_avg = obj_fn(bank)

        with open(f"message.log", "a") as log_f:
            log_f.write(
                f"{round_:>4}   {dcut:4.3f}    {niter:3}   {final_avg.mean():6.3f}   {final_avg.min():6.3f}   "
                f"{final_avg.max():6.3f}    {n_replace:7}   {(time.time() - timechk)/60:8.2f}\n"
            )

        bank[:, 2] = True  # reset to unused solutions

        plot_list.append(final_avg.mean())
        tmp_bank = np.empty([bank.shape[0], 2 + nfeatures], dtype=object)
        tmp_bank[:, 0] = bank[:, 0]
        tmp_bank[:, 1:-1] = bank[:, 3:]
        tmp_bank[:, -1] = final_avg
        tmp_bank[:, 1:] = tmp_bank[:, 1:].astype(np.float16)
        tmp = np.argsort(tmp_bank[:, -1])
        save_bank[round_] = tmp_bank[tmp[::-1]]

    final_bank = np.empty([bank.shape[0], 2 + nfeatures], dtype=object)
    final_bank[:, 0] = bank[:, 0]
    final_bank[:, 1:-1] = bank[:, 3:]
    final_bank[:, -1] = final_avg
    final_bank[:, 1:] = final_bank[:, 1:].astype(np.float16)

    print(f"Total Cost Time: {total_time.get():.3f} min")

    np.save(f"list_bank_{R_d:.5f}_{args.coefficient:.3f}.npy", save_bank)
    save_smiles = pd.DataFrame(save_bank[:, :, 0])
    save_smiles.to_csv(f"list_smiles_{R_d:.5f}_{args.coefficient:.3f}.csv",
                       header=False,
                       index=False)

    df = pd.DataFrame(final_bank, columns=column_name)
    df.to_csv(f"final_bank_{R_d:.5f}_{args.coefficient:.3f}.csv", index=False)

    # log_f.close()
    fail_f.close()

    plt.plot(plot_list)
    plt.tight_layout()
    plt.savefig("target_plot.png")
