import random

explvIn = open('2-explv.bed', 'r')

explv = []
for line in explvIn:
    if line[0] is not '#':
        explv.append(float(line[0:-1].split('\t')[-1]))

#print(explv)

geneIn = open('1-gene.out', 'r')

gene = []
for line in geneIn:
    if 'Gene' not in line and line[0] is not '\n':
        gene.append(line[:-1])

#print(gene)

isoformIn = open('1-transcript.out', 'r')

isoform = []
for line in isoformIn:
    if not 'Gene' in line and line[0] is not '\n':
        substr = line[4:-2].split(', ')
        #print(substr)
        temp = []
        for s in substr:
            temp.append(int(s))
        isoform.append(temp)

#print(isoform)

codes = []
for iso in isoform:
    s = ''
    for id in iso:
        s = s + gene[id]
    codes.append(s)

for i in range(len(codes)):
    print('The ' + str(i) + 'th isoform: ' + str(len(codes[i])))

#print(codes)

#===============================================================================
# N = 100000
# L = 75
# totExp = 0.0
# reads = []
# explv = [1]
# for i in range(len(explv)):
#     totExp = totExp + explv[i] * len(codes[i])
# for i in range(len(explv)):
#     #n = int(N * explv[i] * len(codes[i])/ totExp)
#     n = int(1 * explv[i] * len(codes[i]))
#     print(n)
#     st = 0
#     for j in range(n):
#         #st = random.randint(0, len(codes[i]) - L)
#         read = codes[i][st:st+L]
#         if(len(read) < 75):
#             print('Err!')
#         reads.append(read)
#         st += 1
#         if st + L > len(codes[i]):
#             st = 0
#===============================================================================

reads = []
st = 0
while st + 25 <= len(codes[0]):
    for i in range(10):
        reads.append(codes[0][st:st+25])
    st += 1

outputRead = open('4-naiveReads.fa', 'w')
for x in reads:
    print(x, file = outputRead)
