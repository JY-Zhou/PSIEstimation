from src.KmerHash import KmerHash
from src.EMAlgorithm import EMAlgorithm
import numpy as np 
import json

exonBoundaryFile = r'../input/exonBoundary.bed'
genomeFile = r'../input/genome.fa'
readsFile = r'../input/reads.fq'
K = 15
readLength = 75
kmerHasher = KmerHash(K, readLength, genomeFile, exonBoundaryFile, readsFile)
for x in kmerHasher.temp:
    #print(x)
    print('#'+str(kmerHasher.id[x]), end = '\t')
    if x in kmerHasher.kmerTable:
        print(kmerHasher.kmerTable[x])
    else:
        print('0')
solver = EMAlgorithm(kmerHasher)
solver.work(20)


#===============================================================================
# print('\n======Model Solution==================')
# print(solver.X[0])
# #solver.X[0] = solver.X[0] / (np.ones((1, solver.NX[0])).dot(solver.X[0].T)[0,0])
# print(solver.likelihoodFunction(solver.X[0]))
# print(solver.Psi)
# print('\nConstraints')
# print(solver.A[0].dot(solver.X[0].T))
# print((solver.A[0].dot(solver.X[0].T) > -1e-15).all())
# print(np.ones((1, solver.NX[0])).dot(solver.X[0].T))
# 
# print('\n======Ground Truth====================')
# GX = json.load(open('../kits/XGroundTruth.json', 'r'))
# print(GX)
# print(solver.likelihoodFunction(np.matrix(GX)))
# GPsi = json.load(open('../kits/PsiGroundTruth.json', 'r'))
# print(GPsi)
#===============================================================================
