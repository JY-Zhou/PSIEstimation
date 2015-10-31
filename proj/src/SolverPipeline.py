from src.KmerHash import KmerHash
from src.EMAlgorithm import EMAlgorithm

exonBoundaryFile = r'../input/exonBoundary.bed'
genomeFile = r'../input/genome.fa'
readsFile = r'../input/reads.fq'
K = 25
hasher = KmerHash(K, genomeFile, exonBoundaryFile, readsFile)
solver = EMAlgorithm()
solver.work()