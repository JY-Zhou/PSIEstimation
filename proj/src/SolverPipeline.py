from src.KmerHash import KmerHash
from src.EMAlgorithm import EMAlgorithm

import matplotlib.pyplot as plt
import numpy as np

exonBoundaryFile = r'../input/exonBoundary.bed'
genomeFile = r'../input/genome.fa'
readsFile = r'../input/reads.fq'
K = 25
readLength = 75
kmerHasher = KmerHash(K, readLength, genomeFile, exonBoundaryFile, readsFile)
for x in kmerHasher.temp:
    print(x)
    print('#'+str(kmerHasher.id[x]), end = '\t')
    if x in kmerHasher.kmerTable:
        print(kmerHasher.kmerTable[x])
    else:
        print('0')
solver = EMAlgorithm(kmerHasher)
solver.work(1)

xax = []
y1ax = []
y2ax = []
y3ax = []
vx = np.matrix([[0.4, 0.4, 0.2]])
sx2 = solver.Tau.dot(vx.T)
print(vx)
print(solver.likelihoodFunction(vx))
sx3 = solver.Tau.dot(solver.X[0].T)
print(solver.X[0])
print(solver.likelihoodFunction(solver.X[0]))

i = 0
for x in kmerHasher.temp:
    xax.append(kmerHasher.id[x])
    y1ax.append(kmerHasher.kmerTable[x][0])
    #y2ax.append(sx2[i, 0]*1000)
    #y3ax.append(sx3[i, 0]*1000)
    y2ax.append(y1ax[-1] / sx2[i, 0])
    y3ax.append(y1ax[-1] / sx3[i, 0])
    #y1ax[-1] *= 1000
    i += 1
    
plt.grid(True)
#plt.plot(xax, y1ax, 'b')
plt.plot(xax, y2ax, 'r')
plt.plot(xax, y3ax, 'g')
plt.show()