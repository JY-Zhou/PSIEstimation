import random

class ReferenceGenerator:
    def __init__(self, NG, NE, L, Iso):
        self.NG = NG
        self.NE = NE
        self.L = L
        self.gene = []
        self.genome = []
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

    def generateIntron(self, length):
        ret = ''
        i = 0
        while i < length:
            ret += self.generateBP()
            i += 1
        return ret 

    def generateIntergenic(self, length):
        ret = ''
        i = 0
        while i < length:
            ret += self.generateBP()
            i += 1
        return ret 

    def generateGenome(self):
        self.genome.append(self.generateIntergenic(5000))
        for g in range(self.NG):
            self.gene.append([])
            for e in range(self.NE[g]):
                self.gene[g].append(self.generateExon(self.L[g][e]))
                if not e == self.NE[g]:
                    self.gene[g].append(self.generateIntron(2000))
            self.genome.append(self.generateIntergenic(5000))
        return
       
    def outputRefGenomeFasta(self):
        genomeOut = open('../input/chrSim.fa', 'w')
        s = ''
        for g in range(self.NG):
            s += self.genome[g]
            for i in range(self.NE[g]):
                s += self.gene[g][i*2]
                if not i == self.NE[g]:
                    s += self.gene[g][i*2 + 1]
        s += self.genome[self.NG]

        print('>chrSim', file = genomeOut)
        st = 0
        while len(s[st: st+50]) > 0:
            print(s[st: st+50], file = genomeOut)
            st += 50
    
    def outputExonBoundaryBed(self):
        exonBoundaryOut = open('../input/exonBoundarySim.bed', 'w')
        s = 0
        for g in range(self.NG):
            s += 5000
            st = []
            ed = []
            #for e in range(self.NE[g]):
            #    st.append(s)
            #    s += len(self.gene[g][e])
            #    ed.append(s)
            for i in range(self.NE[g]):
                st.append(s)
                s += len(self.gene[g][i*2])
                ed.append(s)
                if not i == self.NE[g]:
                    s += len(self.gene[g][i*2 + 1])
            print('chrSim', end = '\t', file = exonBoundaryOut)
            print(st[0], end = '\t', file = exonBoundaryOut)
            print(ed[-1], end = '\t', file = exonBoundaryOut)
            print('Gene' + str(g), end = '\t', file = exonBoundaryOut)
            print(0, end = '\t', file = exonBoundaryOut)
            print('+', end = '\t', file = exonBoundaryOut)
            print('.', end = '\t', file = exonBoundaryOut)
            print('.', end = '\t', file = exonBoundaryOut)
            print('0,0,0', end = '\t', file = exonBoundaryOut)

            print(self.NE[g], end = '\t', file = exonBoundaryOut)
            for e in range(self.NE[g]):
                print(ed[e] - st[e], end = ',', file = exonBoundaryOut)
            print('', end = '\t', file = exonBoundaryOut)
            for e in range(self.NE[g]):
                print(st[e] - st[0], end = ',', file = exonBoundaryOut)
            print('', file = exonBoundaryOut)

    def outputExonBoundaryGtf(self):
        exonBoundaryOut = open('../input/exonBoundarySim.gtf', 'w')
        s = 1
        for g in range(self.NG):
            s += 5000
            st = []
            ed = []
            #for e in range(self.NE[g]):
            #    st.append(s)
            #    s += len(self.gene[g][e])
            #    ed.append(s)
            for i in range(self.NE[g]):
                st.append(s)
                s += len(self.gene[g][i*2])
                ed.append(s)
                if not i == self.NE[g]:
                    s += len(self.gene[g][i*2 + 1])
            for e in range(self.NE[g]):
                print('chrSim', end = '\t', file = exonBoundaryOut)
                print('geneGenerator', end = '\t', file = exonBoundaryOut)
                print('exon', end = '\t', file = exonBoundaryOut)
                print(st[e], end = '\t', file = exonBoundaryOut)
                print(ed[e] - 1, end = '\t', file = exonBoundaryOut)
                print(0.0, end = '\t', file = exonBoundaryOut)
                print('+', end = '\t', file = exonBoundaryOut)
                print('.', end = '\t', file = exonBoundaryOut)
                print('gene_id \"Gene' + str(g) + '\"', end = '; ', file = exonBoundaryOut)
                print('transcript_id \"Gene' + str(g) + '\"', end = '; ', file = exonBoundaryOut)
                print('', file = exonBoundaryOut)
            
    def outputIsoformFasta(self):
        isoformOut = open('../kits/isoformSim.fa', 'w')
        for g in range(self.NG):
            for i in range(len(self.Iso[g])):
                isoform = self.Iso[g][i]
                s = ''
                for e in range(len(isoform)):
                    s += self.gene[g][isoform[e]*2]
                print('>Gene' + str(g) + 'Iso' + str(i) + ':W:S', file = isoformOut)
                print(s, file = isoformOut)

    def outputIsoformGtf(self):
        isoformOut = open('../kits/isoformSim.gtf', 'w')
        s = 1
        for g in range(self.NG):
            s += 5000
            st = []
            ed = []
            #for e in range(self.NE[g]):
            #    st.append(s)
            #    s += len(self.gene[g][e])
            #    ed.append(s)
            for i in range(self.NE[g]):
                st.append(s)
                s += len(self.gene[g][i*2])
                ed.append(s)
                if not i == self.NE[g]:
                    s += len(self.gene[g][i*2 + 1])
            k = 0
            for i in range(len(self.Iso[g])):
                iso = self.Iso[g][i]
                for j in range(len(iso)):
                    e = iso[j]
                    print('chrSim', end = '\t', file = isoformOut)
                    print('geneGenerator', end = '\t', file = isoformOut)
                    print('exon', end = '\t', file = isoformOut)
                    print(st[e], end = '\t', file = isoformOut)
                    print(ed[e] - 1, end = '\t', file = isoformOut)
                    print(0.0, end = '\t', file = isoformOut)
                    print('+', end = '\t', file = isoformOut)
                    print('.', end = '\t', file = isoformOut)
                    print('gene_id \"Gene' + str(g) + '\"', end = '; ', file = isoformOut)
                    print('transcript_id \"Gene' + str(g) + 'Iso' + str(k) + '\"', end = '; ', file = isoformOut)
                    print('', file = isoformOut)
                k += 1

    def work(self):
        self.generateGenome()
        self.outputRefGenomeFasta()
        self.outputExonBoundaryBed()
        self.outputExonBoundaryGtf()
        self.outputIsoformFasta()
        self.outputIsoformGtf()
    
