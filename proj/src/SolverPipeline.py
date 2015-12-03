from src.KmerHash import KmerHash
from src.EMAlgorithm import EMAlgorithm
import numpy as np 
import scipy as scp
import json
import time

timeSt = time.clock()

exonBoundaryFile = r'../input/exonBoundary.bed'
genomeFile = r'../input/genome.fa'
readsFile = r'../input/reads.fq'
K = 15
readLength = 75
kmerHasher = KmerHash(K, readLength, genomeFile, exonBoundaryFile, readsFile)
#===============================================================================
# for x in kmerHasher.temp:
#     #print(x)
#     print('#'+str(kmerHasher.id[x]), end = '\t')
#     if x in kmerHasher.kmerTable:
#         print(kmerHasher.kmerTable[x])
#     else:
#         print('0')
#===============================================================================
solver = EMAlgorithm(kmerHasher)
solver.work(20)

print('\n======Model Solution==================')
print('Model\'s X')
print(solver.X)
print('Model\'s Psi')
for val in solver.Psi:
    print(val)
#===============================================================================
# print('\nConstraints')
# print(solver.A[0].dot(solver.X[0].T))
# print((solver.A[0].dot(solver.X[0].T) > -1e-15).all())
# print(np.ones((1, solver.NX[0])).dot(solver.X[0].T))
#===============================================================================
 
print('\n======Ground Truth====================')
GX = json.load(open('../kits/XGroundTruth.json', 'r'))
print('Ground Truth\'s X')
print(GX)
GPsi = json.load(open('../kits/PsiGroundTruth.json', 'r'))
print('Ground Truth\'s Psi')
for val in GPsi:
    print(val)

print('\n======Evaluation============')
print('RMSE:')
mse = 0.0
nmse = 0
for g in range(solver.NG):
    for e in range(solver.NE[g]):
        mse += (solver.Psi[g][e] - GPsi[g][e]) ** 2
        nmse += 1
print(np.sqrt(mse / nmse)) 
print('Time: ' + str(time.clock() - timeSt) + 's')