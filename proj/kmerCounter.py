import os
import sys

#K = int(sys.argv[1])
K = 25
print(K)

filein = open('4-naiveReads.fa', 'r')
#filein = open('4-reads.fa', 'r')
#filein = open('strandedread.fa', 'r')

kmer = {}
cnt = 0
substr = '' 
'''
for line in filein:
    if line[0] is '>':
        if len(substr) is not 0:
            s = substr
            revs = s[::-1]
            substr = ''
            #print(s)
            #print(len(s))
            i = 0
            while i + K <= len(s):
                temp = s[i:i+K]
                if temp in kmer.keys():
                    kmer[temp] = kmer[temp] + 1
                else:
                    kmer[temp] = 1
                #for j in range(i):
                #    print(' ', end='')
                #print(s[i:i+K])
                i = i + 1
    else:
        substr = substr + line[0:-1]

    cnt = cnt + 1
    if cnt % 30000 is 0:
        print(str(cnt//3) + ' lines finished..')
'''
 
totKmer = 0
for line in filein:
    s = line[:-1]
    i = 0
    while i + K <= len(s):
        totKmer += 1
        temp = s[i:i+K]
        if len(temp) < 25:
            print('Err!')
        if temp in kmer.keys():
            kmer[temp] = kmer[temp] + 1
        else:
            kmer[temp] = 1
        i = i + 1
    cnt = cnt + 1
    if cnt % 10000 is 0:
        print(str(cnt) + ' lines finished..')
print('The number of kmers totally is: ' + str(totKmer))

    #input()
#print(len(kmer))
#for x in kmer.keys():
#    print(x + ' ' + str(kmer[x]))

kmerCnt = 0
for x in kmer.keys():
    kmerCnt = kmerCnt + 1
    val = kmer[x]
    kmer[x] = (val, {})
    #print(kmer[x])
print('Kmer number:' + str(kmerCnt))

#input()

gene = []
exonIn = open('1-gene.out', 'r')
exons = []
for line in exonIn:
    #print(line)
    if 'Gene' in line or line is '\n':
        if len(exons) > 0:
            gene.append(exons)
            exons = []
    else:
        exons.append(line[0:-1])

#print('OK!')


miss = 0
total = 0

kmera = []
kmerb = []
for g in range(len(gene)):
    kmera.append([])
    kmerb.append([])
    for i in range(len(gene[g])):
        kmera[g].append(0.0)
        kmerb[g].append([])
        for j in range(len(gene[g])):
            kmerb[g][i].append(0.0)

for g in range(len(gene)):
    for i in range(len(gene[g])):
        label = str(g) + '-' + str(i)
        s = gene[g][i]
        t = 0
        while t + K <= len(s):
            total = total + 1
            temp = s[t:t+K]
            if temp in kmer.keys():
                if label in kmer[temp][1].keys():
                    kmer[temp][1][label] = kmer[temp][1][label] + 1
                else:
                    kmer[temp][1][label] = 1
                kmera[g][i] += kmer[temp][0]
            else:
                print('First type miss!')
                miss = miss + 1
            t = t + 1
        #print(str(t) + '\t' + str(len(s) + 1 - K))

    for i in range(len(gene[g])):
        for j in range(i+1, len(gene[g])):
            #print(str(i) + '\t' + str(j))
            label = str(g) + '-' + str(i) + '-' + str(j)
            s = gene[g][i][-(K-1):] + gene[g][j][0:K-1]
            print(label + ':' + str(len(s)))
            #print(gene[g][i])
            #print(gene[g][j])
            #print(s)
            #print(len(gene[g][i][-(K-1):])) 
            #print(len(gene[g][j][0:K]))
            t = 0
            while t + K <= len(s):
                total = total + 1
                temp = s[t:t+K]
                if temp in kmer.keys():
                    if label in kmer[temp][1].keys():
                        kmer[temp][1][label] = kmer[temp][1][label] + 1
                    else:
                        kmer[temp][1][label] = 1
                    kmerb[g][i][j] += kmer[temp][0]
                else:
                    #print('Second type miss!')
                    miss = miss + 1
                t = t + 1
            #print(t)
print('Unexpressed Frag:' + str(miss) + ', TOTAL: ' + str(total))
#input()

#for x in kmer.keys():
#    print(x)
#    print(kmer[x][0])
#    for y in kmer[x][1].keys():
#        print(y + ': ' + str(kmer[x][1][y]))
#    input()

hashtable = open('5-kmerTable.out', 'w')
z = 0
mul = 0
for key in kmer.keys():
    print(key, end='\t', file=hashtable)
    print(kmer[key][0], end='\t', file=hashtable)
    if len(kmer[key][1].keys()) == 0:
        z = z + 1
        #print(kmer[key][0])
    else:
        num = 0
        for x in kmer[key][1].keys():
            num += kmer[key][1][x]
            print(x, end=',', file=hashtable)
            print(kmer[key][1][x], end=';', file=hashtable)
            if(kmer[key][1][x] != 1):
                print('Not 1 !!')
        if num > 1:
            mul += 1 
            print('num = ' + str(num))
            print('--' + str(kmer[key][1].keys()))
    print('', file=hashtable)
print('Unaligned kmers: ' + str(z))
print('Multi-kmers:' + str(mul))
hashtable.close()

print('---kmer distribution---')
print('kmerA statistic: ')
for g in range(len(kmera)):
    for i in range(len(kmera[g])):
        print(kmera[g][i], end = '\t')
print('\nkmerB statistic: ')
for g in range(len(kmerb)):
    for i in range(len(kmerb[g])):
        for j in range(len(kmerb[g][i])):
            if i < j:
                print(kmerb[g][i][j], end = '\t')
        print()
    
