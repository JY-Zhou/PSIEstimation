import json

class KmerHash:
    def __init__(self, K, readLength, genomeFile, exonBoundaryFile, readsFile):
        self.K = K
        self.readLength = readLength
        self.kmerTable = {}
        self.geneBoundary = []
        self.readReads(readsFile)
        self.readGenome(genomeFile, exonBoundaryFile)
        json.dump(self.kmerTable, open('../output/kmerTable.json', 'w'))
    
    def readReads(self, readsFile):
        fileIn = open(readsFile, 'r')
        proc = 0
        for line in fileIn:
            if proc % 10000 == 0:
                print(str(proc) + ' reads processed...')
            proc += 1
            read = line[:-1]
            st = 0
            while st + self.K <= self.readLength:
                kmer = read[st:st+self.K]
                if kmer in self.kmerTable:
                    self.kmerTable[kmer] += 1
                else:
                    self.kmerTable[kmer] = 1
                st += 1
        self.NW = len(self.kmerTable)
        for kmer in self.kmerTable:
            val = self.kmerTable[kmer]
            self.kmerTable[kmer] = (val, {})
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
        
        self.temp = []
        self.id = {}
        l = 0
        while l + self.K <= len(geneSeq):
            self.temp.append(geneSeq[l:l+self.K])
            self.id[geneSeq[l:l+self.K]] = l
            l += 1
        
        for g in range(self.NG):
            print('Gene-' + str(g) + ' processed...')
            for e in range(self.NE[g]):
                id = str(g) + ',' + str(e) + ',' + str(e)
                st = self.geneBoundary[g][e][0]
                ed = self.geneBoundary[g][e][1] + 1
                
                l = st
                while l + self.K <= ed:
                    kmer = geneSeq[l:l+self.K]
                    if kmer in self.kmerTable:
                        contribution = self.kmerContribution(st, ed, l, l + self.K, ed - st)
                        if id in self.kmerTable[kmer][1]:
                            self.kmerTable[kmer][1][id] += contribution
                        else:
                            self.kmerTable[kmer][1][id] = contribution
                    l += 1
            
                tot = (ed - st - self.readLength + 1) * (self.readLength - self.K + 1)
                    
                l = st
                while l + self.K <= ed:
                    kmer = geneSeq[l:l+self.K]
                    if kmer in self.kmerTable:
                        self.kmerTable[kmer][1][id] /= tot
                    l += 1
                    
            for ei in range(self.NE[g]):
                for ej in range(ei + 1, self.NE[g]):
                    id = str(g) + ',' + str(ei) + ',' + str(ej)
                    edi = self.geneBoundary[g][ei][1] + 1
                    stj = self.geneBoundary[g][ej][0]
                    junction = geneSeq[edi - self.readLength + 1:edi] + geneSeq[stj:stj + self.readLength - 1]
                    
                    l = 0
                    while l + self.K <= 2*self.readLength - 2:
                        kmer = junction[l:l + self.K]
                        if kmer in self.kmerTable:
                            contribution = self.kmerContribution(0, 2*self.readLength - 2, l, l + self.K, 2*self.readLength - 2)
                            if id in self.kmerTable[kmer][1]:
                                self.kmerTable[kmer][1][id] += contribution 
                            else:
                                self.kmerTable[kmer][1][id] = contribution
                        l += 1
                        
                    tot = (2*self.readLength - 2 - self.readLength + 1) * (self.readLength - self.K + 1)
                        
                    l = 0
                    while l + self.K <= 2*self.readLength - 2:
                        kmer = junction[l:l+self.K]
                        if kmer in self.kmerTable:
                            self.kmerTable[kmer][1][id] /= tot
                        l += 1
        return
    
    def kmerContribution(self, st, ed, l, r, L):
        cil = min(L - self.readLength + 1, self.readLength - self.K + 1)
        ret = min(l - st + 1, ed - r + 1)
        #=======================================================================
        # ret = min(ret, cil)
        # ret = min(ret, self.readLength - self.K + 1)
        # ret = min(ret, L - self.readLength + 1)
        # return ret / cil
        #=======================================================================
        return min(ret, cil)
        
    