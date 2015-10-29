import random

class ReferenceGenerator:
    def __init__(self, NG, NE, L, Iso):
        self.NG = NG
        self.NE = NE
        self.L = L
        self.gene = []
        self.Iso = Iso
        
    def generateBP(self):
        return random.choice(['A', 'C', 'G', 'T'])
    
    def generateExon(self, length):
        ret = ''
        i = 0
        while i < length:
            ret += self.generateBP()
            i += 1
        return ret 

    def generateGenome(self):
        for g in range(self.NG):
            self.gene.append([])
            for e in range(self.NE[g]):
                self.gene[g].append(self.generateExon(self.L[g][e]))
        return
       
    def work(self):
        self.generateGenome()
        self.outputRefGenomeFasta()
        self.outputExonBoundaryBed()
        self.outputIsoformFasta()
    
    def outputRefGenomeFasta(self):
        genomeOut = open('../input/genome.fa', 'w')
        s = ''
        for g in range(self.NG):
            for e in range(self.NE[g]):
                s += self.gene[g][e]
        print('chr1>', file = genomeOut)
        st = 0
        while len(s[st: st+50]) > 0:
            print(s[st: st+50], file = genomeOut)
            st += 50
    
    def outputExonBoundaryBed(self):
        exonBoundarydOut = open('../input/exonBoundary.bed', 'w')
        for g in range(self.NG):
            print('Gene' + str(g), end = '\t', file = exonBoundarydOut)
            print(self.NE[g], end = '\t', file = exonBoundarydOut)
            st = []
            ed = []
            s = 0
            for e in range(self.NE[g]):
                st.append(s)
                s += len(self.gene[g][e])
                ed.append(s - 1)
            for e in range(self.NE[g]):
                print(st[e], end = ',', file = exonBoundarydOut)
            print('', end = '\t', file = exonBoundarydOut)
            for e in range(self.NE[g]):
                print(ed[e], end = ',', file = exonBoundarydOut)
            print('', file = exonBoundarydOut)
            
    def outputIsoformFasta(self):
        isoformOut = open('../tricky/isoform.fa', 'w')
        for g in range(self.NG):
            for i in range(len(self.Iso[g])):
                isoform = self.Iso[g][i]
                s = ''
                for e in range(len(isoform)):
                    s += self.gene[g][isoform[e]]
                print(s, file = isoformOut)