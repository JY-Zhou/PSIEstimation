import json
from generator.ReferenceGenerator import ReferenceGenerator
from generator.ReadGenerator import ReadGenerator

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
json.dump(Psi, psiFile)