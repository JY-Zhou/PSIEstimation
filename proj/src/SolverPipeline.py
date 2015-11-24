from src.KmerHash import KmerHash
from src.EMAlgorithm import EMAlgorithm
import json

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

GX = json.load(open('../kits/XGroundTruth.json', 'r'))
print(solver.X)
for g in range(solver.NG):
    tot = 0.0
    for i in range(solver.NX[g]):
        GX[g][i] *= solver.L[g][0,i]
        tot += GX[g][i]
    for i in range(solver.NX[g]):
        GX[g][i] /= tot
        solver.X[g][0, i] = GX[g][i]
print(GX)
print(solver.X)
print(solver.likelihoodFunction(GX[0]))
solver.computePSI()
print(solver.Psi)