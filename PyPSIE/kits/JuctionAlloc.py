import numpy as np

#===============================================================================
# print('!')
# a = np.random.rand(100,10000000)
# print('!!')
# s = a.T
# print('!!!')
#===============================================================================
def showPlot(exon):
    print('Corresponding Reads')
    for x in exon:
        print(x)

    print('Kmer Abundance Distribution')
    k = 25

    dist = []
    for i in range(len(str)):
        dist.append(0)
    l = 0
    r = k
    while r < len(str):
        for x in exon:
            if x[l] != '-' and x[r-1] != '-':
                dist[l] += 1
        l += 1
        r += 1
    for x in dist:
        print(chr(ord('A')+x), end = '')
    print('')
    
    stop = False
    while not stop:
        stop = True
        for i in range(len(dist)):
            if dist[i] == 0:
                stop = False
                print('.', end = '')
            else:
                print(' ', end = '')
            dist[i] -= 1
        print('')
    

str = ''
print(ord('U') - ord('A') + ord('[') - ord('A'))
for i in range(50):
    str += 'a'
for i in range(50):
    str += 'b'
print(str)
step = 75
st = 0
exon1 = []
exon2 = []
junc = []
while len(str[st:st+step]) == step:
    read = ''
    for i in range(st):
        read += '-'
    read += str[st:st+step]
    for i in range(st+step, len(str)):
        read += '-' 
    if 'a' in read and 'b' in read:
        junc.append(read)
    elif 'a' in read:
        exon1.append(read)
    elif 'b' in read:
        exon2.append(read)
    st += 1

showPlot(junc)
showPlot(exon1)
showPlot(exon2)