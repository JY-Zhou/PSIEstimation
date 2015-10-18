import random

def generateBP():
    return random.choice(['A', 'C', 'G', 'T'])

def generateExon(length):
    ret = ''
    for i in range(length):
        ret = ret + generateBP()
    return ret

def generateGene(exonNum):
    exons = []
    for i in range(exonNum):
        #exonLength = random.randint(150,200)
        exonLength = 100 + 20 * i
        exons.append(generateExon(exonLength))
    return exons

def outputGene(gene):
    geneout = open('1-gene.out', 'w')
    for i in range(len(gene)):
        print('Gene' + str(i), file = geneout)
        for exon in gene[i]:
            print(exon, file = geneout)
        print('', file = geneout)

def outputRefFasta(gene):
    fastaout = open('1-reference.fa', 'w')
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

def check(temp, isoform):
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

def generateTranscript(gene, transNum):
    isoform = []
    '''
    if transNum is 1:
        isoform.append(list(range(0, len(exons))))
    else:
        candi = range(1, len(exons)-1)
        for i in range(transNum):
            li = [0, ]
            num = random.randint(1, len(exons)-2)
            temp = random.sample(candi, num)
            for x in temp:
                li.append(x)
            li.append(len(exons)-1)
            li.sort()
            while check(li, isoform):
                li = [0, ]
                num = random.randint(1, len(exons)-2)
                temp = random.sample(candi, num)
                for x in temp:
                    li.append(x)
                li.append(len(exons)-1)
                li.sort()
            isoform.append(li)
    '''
    #isoform.append([0,1,2,3,4,5])
    #isoform.append([0,1,2,4,5])
    isoform.append([0,1,2,3])
    #isoform.append([0,1,3])

    print(isoform)
    return isoform

def outputTranscripts(transcripts):
    tranout = open('1-transcript.out', 'w')
    for i in range(len(transcripts)):
        print('Gene' + str(i), file = tranout)
        for j in range(len(transcripts[i])):
            print(str(j) + '> ', end = '', file = tranout)
            print(transcripts[i][j], file =tranout)

def outputBedFile(gene, transcripts):
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

gene = []
geneNum = 1
for i in range(geneNum):
    #gene.append(generateGene(random.randint(4,8)))
    gene.append(generateGene(4))
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
