class ReadGenerator:
    def __init__(self, expLv, depth, readLength):
        self.expLv = expLv
        self.depth = depth
        self.readLength = readLength
        
    def generateUniformRead(self):
        self.reads = []
        isoformIn = open('../kits/isoform.fa', 'r')
        isoform = []
        for line in isoformIn:
            isoform.append(line[:-1])
        d = 0
        while d < self.depth:
            iso = 0
            for g in range(len(self.expLv)):
                for i in range(len(self.expLv[g])):
                    abund = 0
                    while abund < self.expLv[g][i]:
                        st = 0
                        while len(isoform[iso][st:st+self.readLength]) == self.readLength:
                            self.reads.append(isoform[iso][st:st+self.readLength])
                            st += 1
                        abund += 1
                    iso += 1
            d += 1
        
    def generateProbalisticRead(self):
        pass
    
    def outputReadFastq(self):
        readOut = open('../input/reads.fq', 'w')
        for read in self.reads:
            print(read, file = readOut)
    
    def work(self):
        self.generateUniformRead()
        self.outputReadFastq()