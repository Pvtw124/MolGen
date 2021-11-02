from modSELFIES import get_selfie_chars
from modSELFIES import sanitize_smiles
import selfies
import numpy as np

def avg(lst):
    return sum(lst) / len(lst)

SASscores = []
QEDscores = []
TARGETscores = []
Selfies = []
#Selfies_init = []
lengthSelfies = []
#lengthSelfies_init = []

f = open("final_bank_0.99730_0.900.csv")
#f2 = open("init_bank_0.99730_0.900.csv")
for row in f:
    row = row.split(',')
    Selfies.append(selfies.encoder(row[0]))
    SASscores.append(row[1])
    QEDscores.append(row[2])
    TARGETscores.append(row[3])

# for row in f2:
#     row = row.split(',')
#     Selfies_init.append(selfies.encoder(sanitize_smiles(row[0])))
# print(Selfies_init)
Selfies.pop(0)
#Selfies_init.pop(0)
SASscores.pop(0)
QEDscores.pop(0)
TARGETscores.pop(0)

for i in range(0, len(QEDscores)):
    Selfies[i] = get_selfie_chars(Selfies[i])
    lengthSelfies.append(len(Selfies[i]))
    SASscores[i] = float(SASscores[i])
    QEDscores[i] = float(QEDscores[i])
    TARGETscores[i] = float(TARGETscores[i])
    #lengthSelfies_init.append(len(Selfies_init[i]))




lengthSelfies = np.array(lengthSelfies)
#lengthSelfies_init = np.array(lengthSelfies_init)
    
avgSAS =  avg(SASscores)
avgQED = avg(QEDscores)
avgTARGET = avg(TARGETscores)
print(f" Max length = {np.max(lengthSelfies)}\n Min length = {np.min(lengthSelfies)}\n Average length = {np.mean(lengthSelfies)}\n Standard Deviation = {np.std(lengthSelfies)}\n Average QED score = {avgQED}\n Average TARGET score = {avgTARGET}\n Average SAS score = {avgSAS}")
