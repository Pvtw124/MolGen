def avg(lst):
    return sum(lst) / len(lst)

QEDscores = []
TARGETscores = []

f = open("final_bank_0.99730_0.900.csv")
for row in f:
    row = row.split(',')
    QEDscores.append(row[2])
    TARGETscores.append(row[3])
QEDscores.pop(0)
TARGETscores.pop(0)

for i in range(0, len(QEDscores)):
    QEDscores[i] = float(QEDscores[i])
    TARGETscores[i] = float(TARGETscores[i])
    

avgQED = avg(QEDscores)
avgTARGET = avg(TARGETscores)
print(f"Average QED score = {avgQED}\nAverage TARGET score = {avgTARGET}")
