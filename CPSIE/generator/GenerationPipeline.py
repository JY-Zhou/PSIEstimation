import json
from ReferenceGenerator import ReferenceGenerator
from ReadGenerator import ReadGenerator

confFile = open('../kits/GenerationConfig.json', 'r')

groundTruth = json.load(confFile)
referenceGen = ReferenceGenerator(groundTruth['NG'],
                                  groundTruth['NE'],
                                  groundTruth['L'],
                                  groundTruth['Iso'])
referenceGen.work()
readGen = ReadGenerator(groundTruth['expLv'],
                        groundTruth['depth'],
                        groundTruth['readLength'])
readGen.work()

Psi = []
for g in range(groundTruth['NG']):
    Psi.append([])
    for e in range(groundTruth['NE'][g]):
        Psi[g].append(0.0)
    tot = 0.0
    for i in range(len(groundTruth['Iso'][g])):
        expLv = groundTruth['expLv'][g][i]
        tot += expLv
        for e in range(len(groundTruth['Iso'][g][i])):
            Psi[g][groundTruth['Iso'][g][i][e]] += expLv
    for e in range(groundTruth['NE'][g]):
        Psi[g][e] /= tot

psiFile = open('../kits/PsiGroundTruth.json', 'w')
json.dump(Psi, psiFile, indent = 4)

X = []
for g in range(groundTruth['NG']):
    X.append([])
    dic = {}
    k = 0
    for e in range(groundTruth['NE'][g]):
        X[g].append(0.0)
        k += 1
    for ei in range(groundTruth['NE'][g]):
        for ej in range(ei + 1, groundTruth['NE'][g]):
            X[g].append(0.0)
            dic[(ei, ej)] = k            
            k += 1
            
    for i in range(len(groundTruth['Iso'][g])):
        expLv = groundTruth['expLv'][g][i]
        for e in groundTruth['Iso'][g][i]:
            X[g][e] += expLv
        for a in range(1, len(groundTruth['Iso'][g][i])):
            ei = groundTruth['Iso'][g][i][a-1]
            ej = groundTruth['Iso'][g][i][a]
            X[g][dic[(ei, ej)]] += expLv
            
    tot = 0.0
    for x in X[g]:
        tot += x
    for x in range(len(X[g])):
        X[g][x] /= tot
XFile = open('../kits/XGroundTruth.json', 'w')
json.dump(X, XFile, indent = 4)