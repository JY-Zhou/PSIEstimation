import random

class ReferGenerator:
    def __init__(self, NG, NE, L):
        self.NG = NG
        self.NE = NE
        self.L = L        
        self.gene = []
        
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
            for e in range(self.NE):
                self.gene[g].append(self.generateExon(self.L[g][e]))
        return
       
    def generateTranscript(self, ):
        isoform = []
        isoform.append([0,1,2,3])
        print(isoform)
        return isoform

    def outputRefFasta(self, gene):
        fastaout = open('../input/reference.fa', 'w')
        s = ''
        for i in range(len(gene)):
            for exon in gene[i]:
                s = s + exon
        print('>chr1', file = fastaout)
        k = 0
        for ch in s:
            print(ch, end = '', file = fastaout)
            k = k + 1
            if k == 50:
                print('', file = fastaout)
                k = 0
                
    def outputGene(self, gene):
        geneout = open('1-gene.out', 'w')
        for i in range(len(gene)):
            print('Gene' + str(i), file = geneout)
            for exon in gene[i]:
                print(exon, file = geneout)
            print('', file = geneout)


    def check(self, temp, isoform):
        for x in isoform:
            if len(x) == len(temp):
                flag = True
                for i in range(len(x)):
                    if x[i] != temp[i]:
                        flag = False
                        break
                if flag:
                    return True
        return False

    

    def outputTranscripts(self, transcripts):
        tranout = open('1-transcript.out', 'w')
        for i in range(len(transcripts)):
            print('Gene' + str(i), file = tranout)
            for j in range(len(transcripts[i])):
                print(str(j) + '> ', end = '', file = tranout)
                print(transcripts[i][j], file =tranout)

    def outputBedFile(self, gene, transcripts):
        outfile = open('1-sample.bed', 'w')
        exonsite = open('1-exonlocus.out', 'w')
        genest = []
        geneed = []
        k = 0
        for i in range(len(gene)):
            exonst = []
            exoned = []
            exonst.append(k)
            for j in range(len(gene[i])):
                k = k + len(gene[i][j])
                exoned.append(k)
                if j < len(gene[i]) - 1:
                    exonst.append(k)
            print('Gene' + str(i), file = exonsite)
            print(exonst, file = exonsite)
            print(exoned, file = exonsite)
            genest.append(exonst)
            geneed.append(exoned)

        for i in range(len(transcripts)):
            gst = genest[i][0]
            ged = geneed[i][len(geneed[i]) - 1]
            for j in range(len(transcripts[i])):
                name = 'gene' + str(i) + '.' + 'isoform' + str(j)
                score = '0'
                strand = '+'
                color = '0'
                blockNum = str(len(transcripts[i][j]))
                blockSize = ''
                blockSt = ''
                for p in range(len(transcripts[i][j])):
                    exonid = transcripts[i][j][p]
                    blockSize = blockSize + str(len(gene[i][exonid])) + ','
                    blockSt = blockSt + str(genest[i][exonid] - genest[i][0]) + ','
                line = 'chr1\t' + str(gst) + '\t' + str(ged) + '\t';
                line = line + name + '\t' + score + '\t' + strand + '\t'
                line = line + str(gst) + '\t' + str(ged) + '\t';
                line = line + color + '\t' + blockNum + '\t'
                line = line + blockSize + '\t' + blockSt
                print(line, file = outfile)

if __name__ == '__main__':
    gen = ReferGenerator()
    gene = []
    geneNum = 1
    for i in range(geneNum):
        #gene.append(generateGene(random.randint(4,8)))
        gene.append(gen.generateGene(4))
        #gene.append(generateGene(6))

    transcripts = []
    for i in range(len(gene)):
        #transcripts.append(generateTranscript(gene[i], random.randint(1,3)))
        #transcripts.append(generateTranscript(gene[i], 2))
        transcripts.append(generateTranscript(gene[i], 1))

    outputGene(gene)
    outputRefFasta(gene)
    outputTranscripts(transcripts)
    outputBedFile(gene, transcripts)

    for i in range(len(gene)):
        print('Gene' + str(i) + ': ' + str(len(transcripts[i])))
