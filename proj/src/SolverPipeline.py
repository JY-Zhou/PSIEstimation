from src.KmerHash import KmerHash
from src.EMAlgorithm import EMAlgorithm

exonBoundaryFile = r'../input/exonBoundary.bed'
genomeFile = r'../input/genome.fa'
readsFile = r'../input/reads.fq'
K = 25
readLength = 75
kmerHasher = KmerHash(K, readLength, genomeFile, exonBoundaryFile, readsFile)
for x in kmerHasher.temp:
    print(x)
    print('#'+str(kmerHasher.contri[x]), end = '\t')
    print(kmerHasher.kmerTable[x])
solver = EMAlgorithm(kmerHasher)
solver.work(1)