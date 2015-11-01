from src.KmerHash import KmerHash
from src.EMAlgorithm import EMAlgorithm

exonBoundaryFile = r'../input/exonBoundary.bed'
genomeFile = r'../input/genome.fa'
readsFile = r'../input/reads.fq'
K = 25
readLength = 75
kmerHasher = KmerHash(K, readLength, genomeFile, exonBoundaryFile, readsFile)
solver = EMAlgorithm()
solver.work()