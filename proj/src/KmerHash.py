from src.SolverPipeline import exonBoundaryFile
class KmerHash:
    def __init__(self, K, genomeFile, exonBoundaryFile, readsFile):
        self.K = K
        self.kmerTable = {}
        self.gene = []
        self.readReads(readsFile)
        self.readGenome(genomeFile, exonBoundaryFile)
    
    def readReads(self, readsFile):
        fileIn = open(readsFile, 'r')
        for line in fileIn:
            read = line[:-1]
            if read in self.kmerTable.keys():
                self.kmerTable[read] += 1
            else:
                self.kmerTable[read] = 0
        self.NW = len(self.kmerTable)
        return 
    
    def readGenome(self, genomeFile, exonBoundaryFile):
        exonBoundaryIn = open(exonBoundaryFile,'r')
        COL_EXONNUM = 1
        COL_EXONST = 2
        COL_EXONED = 3
                
        for line in exonBoundaryIn:
            sub = line[:-1].split('\t')
            self.gene.append([])
            subst = sub[COL_EXONST].split(',')
            subed = sub[COL_EXONED].split(',')
            e = 0
            while e < int(sub[COL_EXONNUM]):
                self.gene[-1].append((int(subst[e]), int(subed[e])))
                e += 1
        
        self.NG = len(self.gene)
        self.NE = []
        for g in range(self.NG):
            self.NE.append(len(self.gene[g]))
        
        exonBoundaryIn.close()
                
        genomeIn = open(genomeFile, 'r')

        #go on ...
                
        
        
        return
    
        
    def work(self):
        return       
        