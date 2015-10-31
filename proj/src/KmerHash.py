class KmerHash:
    def __init__(self, K, genomeFile, exonBoundaryFile, readsFile):
        self.K = K
        self.kmerTable = {}
        self.geneBoundary = []
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
            self.geneBoundary.append([])
            subst = sub[COL_EXONST].split(',')
            subed = sub[COL_EXONED].split(',')
            e = 0
            while e < int(sub[COL_EXONNUM]):
                self.geneBoundary[-1].append((int(subst[e]), int(subed[e])))
                e += 1
        exonBoundaryIn.close()
        
        self.NG = len(self.geneBoundary)
        self.NE = []
        for g in range(self.NG):
            self.NE.append(len(self.geneBoundary[g]))
        
        genomeIn = open(genomeFile, 'r')
        geneSeq = ''
        for line in genomeIn:
            if '>' not in line:
                geneSeq += line[:-1]
        genomeIn.close()
        for g in range(self.NG):
            
            ...
            
        return