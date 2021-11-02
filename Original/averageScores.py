import selfies
import numpy as np

def avg(lst):
    return sum(lst) / len(lst)

def get_selfie_chars(selfie):
    chars_selfie = [] # A list of all SELFIE sybols from string selfie
    while selfie != '':
        chars_selfie.append(selfie[selfie.find('['): selfie.find(']')+1])
        selfie = selfie[selfie.find(']')+1:]
    return chars_selfie

SASscores = []
QEDscores = []
TARGETscores = []
Selfies = []
lengthSelfies = []

f = open("final_bank_0.99730_0.900.csv")
for row in f:
    row = row.split(',')
    Selfies.append(selfies.encoder(row[0]))
    SASscores.append(row[1])
    QEDscores.append(row[2])
    TARGETscores.append(row[3])
Selfies.pop(0)
SASscores.pop(0)
QEDscores.pop(0)
TARGETscores.pop(0)

for i in range(0, len(QEDscores)):
    Selfies[i] = get_selfie_chars(Selfies[i])
    lengthSelfies.append(len(Selfies[i]))
    SASscores[i] = float(SASscores[i])
    QEDscores[i] = float(QEDscores[i])
    TARGETscores[i] = float(TARGETscores[i])


lengthSelfies = np.array(lengthSelfies)
    
avgSAS =  avg(SASscores)
avgQED = avg(QEDscores)
avgTARGET = avg(TARGETscores)
print(f" Max length = {np.max(lengthSelfies)}\n Min length = {np.min(lengthSelfies)}\n Average length = {np.mean(lengthSelfies)}\n Standard Deviation = {np.std(lengthSelfies)}\n Average QED score = {avgQED}\n Average TARGET score = {avgTARGET}\n Average SAS score = {avgSAS}")

