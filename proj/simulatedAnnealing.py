import random
import math
import logging
import EnhancedSA

checkTime = 0

def Q(length, kmer, alpha, beta, K):
    ret = 0.0
    test = 0.0
    bucketa = []
    bucketb = []
    kmera = []
    kmerb = []
    checkera = []
    checkerb = []
     
    for g in range(len(alpha)):
        bucketa.append([])
        bucketb.append([])
        kmera.append([])
        kmerb.append([])
        checkera.append([])
        checkerb.append([])
        for i in range(len(alpha[g])):
            bucketa[g].append(0)
            bucketb[g].append([])
            kmera[g].append(0)
            kmerb[g].append([])
            checkera[g].append(0)
            checkerb[g].append([])
            for j in range(len(beta[g][i])):
                bucketb[g][i].append(0)
                kmerb[g][i].append(0)
                checkerb[g][i].append(0)

    for k in kmer.keys():
        w = kmer[k][0]
        s = 0.0
        newS = []
        newL = []

        #print('kmer info ' + str(len(kmer[k][1].keys())))
        for x in kmer[k][1].keys():
            substr = x.split('-')
            if len(substr) == 2:
                g = int(substr[0])
                i = int(substr[1])
                temp =  kmer[k][1][x] * alpha[g][i] / length[g][i]
                s = s + temp
                newS.append(temp)
                newL.append(length[g][i])
                bucketa[g][i] += 1
                kmera[g][i] += w
                #print('tau: ' + str(kmer[k][1][x]/length[g][i]) + ' occurence: ' + str(kmer[k][1][x]) + ' length~: ' + str(length[g][i]))
                
            elif len(substr) == 3:
                g = int(substr[0])
                i = int(substr[1])
                j = int(substr[2])
                temp = kmer[k][1][x] * beta[g][i][j] / (K - 1)
                s = s + temp
                newS.append(temp)
                newL.append(K-1)
                bucketb[g][i][j] += 1
                kmerb[g][i][j] += w
                #print('tau: ' + str(kmer[k][1][x]/length[g][i]) + ' occurence: ' + str(kmer[k][1][x]) + ' length~: ' + str(K-1))
        #print('s = ' + str(s))
        #if(s <= 0.0):
            #continue
        test += s
        if s > 0:
            s = math.log(s)
        else:
            s = -float('inf')
        #err = 0.0
        #for i in range(len(newL)):
        #    if newS[i] > 0:
        #        err += 1.0 / newL[i] 
        #        #print("Hit!" + str(k))
        #s = w * s * err
        s *= w
        ret += s
    #print(test)
    #print(ret)

  #=============================================================================
  #   print('BucketA statistic: ')
  #   for g in range(len(bucketa)):
  #       for i in range(len(bucketa[g])):
  #           print(bucketa[g][i], end = '\t')
  #   print('\nBucketB statistic: ')
  #   for g in range(len(bucketb)):
  #       for i in range(len(bucketb[g])):
  #           for j in range(len(bucketb[g][i])):
  #               if i < j:
  #                   print(bucketb[g][i][j], end = '\t')
  #           print()
  #     
  #   print('kmerA statistic: ')
  #   for g in range(len(kmera)):
  #       for i in range(len(kmera[g])):
  #           print(kmera[g][i], end = '\t')
  #   print('\nkmerB statistic: ')
  #   for g in range(len(kmerb)):
  #       for i in range(len(kmerb[g])):
  #           for j in range(len(kmerb[g][i])):
  #               if i < j:
  #                   print(kmerb[g][i][j], end = '\t')
  #           print()
  # 
  #   print('normalized kmerA statistic: ')
  #   for g in range(len(kmera)):
  #       for i in range(len(kmera[g])):
  #           print(kmera[g][i] / bucketa[g][i], end = '\t')
  #   print('\nnormalized kmerB statistic: ')
  #   for g in range(len(kmerb)):
  #       for i in range(len(kmerb[g])):
  #           for j in range(len(kmerb[g][i])):
  #               if i < j:
  #                   if bucketb[g][i][j] > 0:
  #                       print(kmerb[g][i][j] / bucketb[g][i][j], end = '\t')
  #                   else:
  #                       print('0', end = '\t')
  #           print()
  # 
  #   resChecker = 0.0
  #   for g in range(len(alpha)):
  #       for i in range(len(alpha[g])):
  #           if alpha[g][i] > 0:
  #               checkera[g][i] = kmera[g][i] * math.log(alpha[g][i] / length[g][i])
  #           else:
  #               checkera[g][i] = -float('inf')
  #           resChecker += checkera[g][i]
  #           for j in range(len(beta[g][i])):
  #               if i < j:
  #                   if beta[g][i][j] > 0:
  #                       checkerb[g][i][j] = kmerb[g][i][j] * math.log(beta[g][i][j] / (K - 1))
  #                   else:
  #                       checkerb[g][i][j] = -float('inf')
  #                   resChecker += checkerb[g][i][j]
  # 
  #   print('Check result = ' + str(resChecker))
  #     
  #   print('checker statistic: ')
  #   for g in range(len(checkera)):
  #       for i in range(len(checkera[g])):
  #           print(checkera[g][i], end = '\t')
  #   print('\nchecker statistic: ')
  #   for g in range(len(checkerb)):
  #       for i in range(len(checkerb[g])):
  #           for j in range(len(checkerb[g][i])):
  #               if i < j:
  #                   print(checkerb[g][i][j], end = '\t')
  #           print()
  #=============================================================================
    return ret

def initial():
    alpha = []
    beta = []
    for g in length:
        alpha.append([])
        beta.append([])
        for x in range(len(g)):
            alpha[-1].append(1.0)
            #alpha[-1].append(random.uniform(0, 1))
            beta[-1].append([])
            for y in range(len(g)):
                if x >= y:
                    beta[-1][x].append(0.0)
                else:
                    beta[-1][x].append(0.1)
                    #beta[-1][x].append(random.uniform(0, 1))

    alpha[-1] = [76, 96, 116, 136]
    beta[-1] = [[0, 20, 0, 0], [0, 0, 24, 0], [0, 0, 0, 24], [0, 0, 0, 0]]
    #alpha[-1] = [100000, 100000, 100000, 100000]
    #beta[-1] = [[0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0]]
    #alpha[-1] = [133, 149, 44, 121]
    #beta[-1] = [[0, 27, 0, 0], [0, 0, 6, 21], [0, 0, 0, 6], [0, 0, 0, 0]]
    #alpha[-1] = [142, 131, 172, 133]
    #beta[-1] = [[0, 24, 0, 0], [0, 0, 24, 24], [0, 0, 0, 24], [0, 0, 0, 0]]
    return normalize(alpha, beta)
    #return (alpha, beta)

def check(alpha, beta, length, K=25):
    global checkTime
    checkTime = checkTime + 1
    showTrend(alpha, beta)
    s = 0
    for g in alpha:
        for i in g:
            s = s + i
    for g in beta:
        for i in range(len(g)):
            for j in range(i+1, len(g[i])):
                if i < j:
                    s = s + g[i][j]
    if abs(s - 1.0) > 1e-8:
        logging.info('Check: Unnomalized')
        return False
    
    for g in range(len(alpha)):
        sa = 0.0
        sb = 0.0
        for i in range(len(alpha[g])):
            sa = sa + alpha[g][i]/length[g][i]
        for i in range(len(beta[g])):
            for j in range(len(beta[g][i])):
                if i < j:
                    sb = sb + beta[g][i][j]/(K-1)
        for i in range(len(alpha[g])):
            if alpha[g][i]/length[g][i] > sa - sb:
                logging.info('Check: Psi over 1')
                return False
        for i in range(len(beta[g])):
            l = 0.0
            r = 0.0
            for j in range(len(beta[g][i])):
                if i > j:
                    l = l + beta[g][j][i]/(K-1)
                if i < j:
                    r = r + beta[g][i][j]/(K-1)
            if alpha[g][i] / length[g][i] < l:
                logging.info('Check: left over')
                return False
            if alpha[g][i] / length[g][i] < r:
                logging.info('Check: right over')
                return False
    return True

def normalize(al, be):
    s = 0
    alpha = []
    beta = []
    for g in range(len(al)):
        alpha.append([])
        for i in range(len(al[g])):
            alpha[g].append(al[g][i])
            s = s + al[g][i]
    for g in range(len(be)):
        beta.append([])
        for i in range(len(be[g])):
            beta[g].append([])
            for j in range(len(be[g][i])):
                beta[g][i].append(be[g][i][j])
                s = s + be[g][i][j]
    for g in range(len(alpha)):
        for i in range(len(alpha[g])):
            alpha[g][i] = alpha[g][i] / s
    for g in range(len(beta)):
        for i in range(len(beta[g])):
            for j in range(len(beta[g][i])):
                beta[g][i][j]  = beta[g][i][j] / s
    return (alpha, beta)
    #return (al, be)

def interfer(al, be, low, up):
    alpha = []
    beta = []
    for g in range(len(al)):
        alpha.append([])
        for i in range(len(al[g])):
            alpha[g].append(al[g][i] * random.uniform(low, up))
    for g in range(len(be)):
        beta.append([])
        for i in range(len(be[g])):
            beta[g].append([])
            for j in range(len(be[g][i])):
                beta[g][i].append(be[g][i][j] * random.uniform(low, up))
 
    return normalize(alpha, beta)

def showTrend(alpha, beta):
    for g in range(len(alpha)):
        for i in range(len(alpha[g])):
            print(alpha[g][i], end = '\t', file = fig)
    for g in range(len(beta)):
        for i in range(len(beta[g])):
            for j in range(len(beta[g][i])):
                if i < j:
                    print(beta[g][i][j], end = '\t', file = fig)
    print('', file = fig)

def show(alpha, beta):
    print('-------Alpha-------')
    for g in alpha:
        for i in g:
            print(i, end='\t')
    print('\n-------Beta--------')
    for g in beta:
        for i in range(len(g)):
            for j in range(len(g[i])):      
                if i < j:
                    print(g[i][j], end='\t')
            print('')
    print('Press any key to continue..')
    #input()

def annealing(length, kmer, K):
    global checkTime
    checkTime = 0
    (alpha, beta) = initial()
    while not check(alpha, beta, length):
        (alpha, beta) = initial()

    print('Check ' + str(checkTime) + ' times to initialize' )
    show(alpha, beta)

    T = 100
    tau = 0.99
    L = 100 
    step = 0.1
    low = 1 - step
    up = 1 + step
    while T > 20:
        l = 0
        f = Q(length, kmer, alpha, beta, K)
        print(str(T) + ': ' + str(f))
        #return (alpha, beta)
        #exit(0)

        cnt = 0
        while l < L:
            cnt += 1
            if cnt % 1000000 == 0:
                print('Loop for :' + str(cnt//1000000) + 'M times')
            (newalpha, newbeta) = interfer(alpha, beta, low, up)
            while not check(newalpha, newbeta, length):
                (newalpha, newbeta) = interfer(alpha, beta, low, up)
            newf = Q(length, kmer, newalpha, newbeta, K)
            #input()
            delta = newf - f
            if delta > 0:
                (alpha, beta) = (newalpha, newbeta)
                f = newf
                l = 0
            else:
                #p = math.exp(delta/(T))
                p = 0 
                if random.uniform(0.0, 1.0) < p:
                    (alpha, beta) = (newalpha, newbeta)
                    f = newf
            l += 1
        T = T * tau
    return (alpha, beta)

def calpsi(alpha, beta, length, k):
    psi = []
    for g in range(len(alpha)):
        psi.append([])
        sa = 0.0
        sb = 0.0
        print('--------Alpha-------')
        for i in range(len(alpha[g])):
            print(alpha[g][i], end = '\t')
            sa = sa + alpha[g][i] / length[g][i]
            #sa = sa + alpha[g][i]
            

        print('\n--------Beta-------')
        for i in range(len(beta[g])):
            for j in range(len(beta[g][i])):
                if i < j:
                    print(beta[g][i][j], end = '\t')
                    sb = sb + beta[g][i][j]/(K-1)
                    #sb = sb + beta[g][i][j]
            print('') 
        
        print('Sum of alpha: ' + str(sa))
        print('Sum of beta:  ' + str(sb))
        for i in range(len(alpha[g])):
            psi[g].append(alpha[g][i] / length[g][i] / (sa - sb))
    return psi

fig = open('ABTrend.out', 'w')
inFile = open('5-kmerTable.out', 'r')

kmer = {}
K = 0
for line in inFile:
    substr = line[:-1].split('\t')
    K = len(substr[0])
    kmer[substr[0]] = (float(substr[1]), {})
    dic = substr[-1][:-1].split(';')
    if len(dic) == 2:
        print('Multi-aligned kmer')
    for x in dic:
        kv = x.split(',')
        kmer[substr[0]][1][kv[0]] = float(kv[1])

length = []
exonIn = open('1-gene.out', 'r')
exons = []
for line in exonIn:
    if 'Gene' in line or line is '\n':
        if len(exons) > 0:
            length.append(exons)
            exons = []
    else:
        exons.append(len(line[:-1]) - K + 1)

(alpha, beta) = annealing(length, kmer, K)

psi = calpsi(alpha, beta, length, K)
output = open('6-predict_SA.psi', 'w')
for g in range(len(psi)):
    print('Gene' + str(g), file=output)
    print('Gene' + str(g))
    print(psi[g], file=output)
    print(psi[g])
