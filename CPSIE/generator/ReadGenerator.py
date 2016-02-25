class ReadGenerator:
    def __init__(self, expLv, depth, readLength):
        self.expLv = expLv
        self.depth = depth
        self.readLength = readLength
        
    def generateUniformRead(self):
        self.reads = []
        isoformIn = open('../kits/isoformSim.fa', 'r')
        isoform = []
        k = 0
        for line in isoformIn:
            if k % 2 == 1:
                isoform.append(line[:-1])
            k += 1

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
        readOut = open('../input/readsSim.fq', 'w')
        k = 0
        for read in self.reads:
            print('@Read' + str(k) +':W:S', file = readOut)
            print(read, file = readOut)
            print('+', file = readOut)
            for i in range(len(read)):
                print('~', end = '', file = readOut)
            print('', file = readOut)
            k += 1
    
    def outputReadFasta(self):
        readOut = open('../input/readsSim.fa', 'w')
        k = 0
        for read in self.reads:
            print('>Read' + str(k) +':W:S', file = readOut)
            print(read, file = readOut)
            k += 1
        

    def work(self):
        self.generateUniformRead()
        self.outputReadFastq()
        self.outputReadFasta()
