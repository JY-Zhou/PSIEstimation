import sys
import math
from EnhancedSA import EnhancedSA
import random

class OptimizeMLE(EnhancedSA):
    RoStepVal = 0.01
    plainStay = 50
    
    # Specific parameters
    K = 0
    kmer = {}
    length = []
    mergeIndex = []
    splitIndex = []
    
    wVec = []
    tauMat = []
        
    def initializeSpecificParameters(self):
        #Read kmer table
        kmerTable = open('5-kmerTable.out', 'r')
        for line in kmerTable:
            substr = line[:-1].split('\t')
            self.K = len(substr[0])
            self.kmer[substr[0]] = (float(substr[1]), {})
            dic = substr[-1][:-1].split(';')
            if len(dic) is 2:
                print('Alert: Includes multi-kmer', file = sys.stderr)
            for entry in dic:
                key = entry.split(',')[0]
                value = entry.split(',')[1]
                self.kmer[substr[0]][1][key] = float(value)
        
        #Read exon lengths
        exonTable = open('1-gene.out', 'r')
        exons = []
        for line in exonTable:
            if 'Gene' in line or line is '\n':
                if len(exons) > 0:
                    self.length.append(exons)
                    exons = []
            else:
                exons.append(len(line[:-1]))
        
        #Calculate the index map
        for g in range(len(self.length)):
            self.mergeIndex.append([])
            for i in range(len(self.length[g])):
                self.mergeIndex[g].append([])
                for j in range(len(self.length[g])):
                    self.mergeIndex[g][i].append([])
        self.n = 0
        for g in range(len(self.mergeIndex)):
            for i in range(len(self.mergeIndex[g])):
                self.mergeIndex[g][i][i] = self.n
                self.splitIndex.append((g, i, i))
                self.n += 1
            for i in range(len(self.mergeIndex[g])):
                for j in range(len(self.mergeIndex[g][i])):
                    if i < j:
                        self.mergeIndex[g][i][j] = self.n
                        self.splitIndex.append((g, i, j))
                        self.n += 1
        print('Check: Totally ' + str(self.n) + ' variables', file = sys.stderr)
        
        #Initial coeffients
        for key in self.kmer.keys():
            self.wVec.append(self.kmer[key][0])
            self.tauMat.append([])
            for i in range(self.n):
                self.tauMat[-1].append(0.0)
            for x in self.kmer[key][1].keys():
                substr = x.split('-')
                if len(substr) is 2:
                    g = int(substr[0])
                    i = int(substr[1])
                    val = self.kmer[key][1][x] / (self.length[g][i] - self.K + 1)
                    self.tauMat[-1][self.mergeIndex[g][i][i]] = val 
                else:
                    g = int(substr[0])
                    i = int(substr[1])
                    j = int(substr[2])
                    val = self.kmer[key][1][x] / (self.K - 1)
                    self.tauMat[-1][self.mergeIndex[g][i][j]] = val
            
        #=======================================================================
        # print(self.wVec)
        # for i in range(len(self.tauMat)):
        #     print(self.tauMat[i])
        #=======================================================================
        
    def boundaryChecker(self, x):
        s = 0
        for v in x:
            s += v
        if math.fabs(s - 1.0) > self.epsRel:
            print('Alert: Failed on checking boundaries (not normalized)', file = sys.stderr)
            return False 
        for g in range(len(self.length)):
            sa = 0.0
            sb = 0.0
            for i in range(len(self.length[g])):
                sa += x[self.mergeIndex[g][i][i]] / (self.length[g][i] - self.K + 1)
                for j in range(len(self.length[g])):
                    if i < j:
                        sb += x[self.mergeIndex[g][i][j]] / (self.K - 1)
            for i in range(len(self.length[g])):
                if x[self.mergeIndex[g][i][i]] / (self.length[g][i] - self.K + 1) > sa - sb:
                    print('Alert: Failed on checking boundaries (psi will higher than 1)', file = sys.stderr)
                    return False
        
            for i in range(len(self.length[g])):
                sl = 0.0
                sr = 0.0
                for j in range(len(self.length[g])):
                    if j < i:
                        sl += x[self.mergeIndex[g][j][i]] / (self.K - 1)
                    if i < j:
                        sr += x[self.mergeIndex[g][i][j]] / (self.K - 1)
                if x[self.mergeIndex[g][i][i]] / (self.length[g][i] - self.K + 1) < sr:
                    print('Alert: Failed on checking boundaries (the abundance of downstream junction invalid)', file = sys.stderr)
                    return False
                if x[self.mergeIndex[g][i][i]] / (self.length[g][i] - self.K + 1) < sl:
                    print('Alert: Failed on checking boundaries (the abundance of upstream junction invalid)', file = sys.stderr)
                    print('--' + str(g) + '--' + str(i) + '< ' + str(sl))
                    return False
        return True       
        
    def normalize(self, x):
        s = 0.0
        for v in x:
            s += v
        for i in range(len(x)):
            x[i] /= s
        return x
            
    def f(self, x):
        ret = 0.0
        for key in self.kmer.keys():
            w = self.kmer[key][0]
            s = 0.0
            for info in self.kmer[key][1].keys():
                y = self.kmer[key][1][info]
                substr = info.split('-')
                if len(substr) is 2:
                    g = int(substr[0])
                    i = int(substr[1])
                    s += y * x[self.mergeIndex[g][i][i]] / (self.length[g][i] - self.K + 1)
                elif len(substr) is 3:
                    g = int(substr[0])
                    i = int(substr[1])
                    j = int(substr[2])
                    s += y * x[self.mergeIndex[g][i][j]] / (self.K - 1)
            if s > 0:
                s = math.log(s)
            else:
                s = -float('inf')
            ret += w * s
        return ret
    
    def getInitialX(self):
        for g in range(len(self.length)):
            for i in range(len(self.length[g])):
                self.x.append(1)
            for i in range(len(self.length[g])):
                for j in range(len(self.length[g])):
                    if i < j:
                        self.x.append(self.epsRel)  
        self.x = self.normalize(self.x)
        print('Check: initial x\'s = ' + str(self.x), file = sys.stderr)
        
        for i in range(self.n):
            self.xMax.append(1.0)
            self.xMin.append(0.0)
            
        if self.boundaryChecker(self.x):
            self.updateBoundary()
        else:
            print('Error: initial failed')
    
        print('Check: initial xMin\'s = ' + str(self.xMin), file = sys.stderr)
        print('Check: initial xMax\'s = ' + str(self.xMax), file = sys.stderr)
        
    def updateBoundary(self):
        scale = []
        for i in range(self.n):
            scale.append(self.xMax[i] - self.xMin[i])

        for g in range(len(self.length)):
            sa = 0.0
            sb = 0.0
            for i in range(len(self.length[g])):
                sa += self.x[self.mergeIndex[g][i][i]] / (self.length[g][i] - self.K + 1)
                for j in range(len(self.length[g])):
                    if i < j:
                        sb += self.x[self.mergeIndex[g][i][j]] / (self.K - 1)
            for i in range(len(self.length[g])):
                self.xMax[self.mergeIndex[g][i][i]] = (sa - sb) * (self.length[g][i] - self.K + 1)
        
            for i in range(len(self.length[g])):
                sl = 0.0
                sr = 0.0
                for j in range(len(self.length[g])):
                    if j < i:
                        sl += self.x[self.mergeIndex[g][j][i]] / (self.K - 1)
                    if i < j:
                        sr += self.x[self.mergeIndex[g][i][j]] / (self.K - 1)
                self.xMin[self.mergeIndex[g][i][i]] = max(sl, sr) * (self.length[g][i] - self.K + 1) 
                
            for i in range(len(self.length[g])):
                for j in range(len(self.length[g])):
                    if i < j:
                        up = min(self.x[self.mergeIndex[g][i][i]], self.x[self.mergeIndex[g][j][j]])
                        self.xMax[self.mergeIndex[g][i][j]] = up * (self.K - 1)                     
                        self.xMin[self.mergeIndex[g][i][j]] = 0.0
                        
        if len(self.step) is self.n:
            for i in range(self.n):
                self.step[i] *= (self.xMax[i] - self.xMin[i]) / scale[i]
                
        for i in range(self.n):
            if self.x[i] < self.xMin[i]:
                self.x[i] = self.xMin[i] + self.epsRel
            if self.x[i] > self.xMax[i]:
                self.x[i] = self.xMax[i] - self.epsRel
    
    def move(self):
        xTry = []
        self.spacePartition()
        for i in range(self.n):
            xInit = self.x[i]
            rand = random.uniform(-1.0, 1.0)
            xNew = xInit
            if i in self.moveIndex:
                xNew = xInit + rand * self.step[i]
                if xNew > self.xMax[i] or xNew < self.xMin[i]:
                    xNew = xInit - rand * self.step[i]
            xTry.append(xNew)
            if xTry[i] > self.xMax[i] or xTry[i] < self.xMin[i]:
                print('Error: Move out of boundaries.' , file = sys.stderr)
                print('--xMin = ' + str(self.xMin[i]) + ', xMax = ' + str(self.xMax[i]) +', index = ' + str(self.splitIndex[i]), file = sys.stderr)
                print('--ValInit = ' + str(xInit) + ' , ValNew = ' + str(xNew) + ', Step = ' + str(rand * self.step[i]), file = sys.stderr)
                print('--Rand = ' + str(rand) + ' , perStep = ' + str(self.step[i]), file = sys.stderr)
                input()
        
        if len(xTry) != self.n:
            print('Error: Dimension loss when moving.', file = sys.stderr)
        
        xTry = self.normalize(xTry)
        if self.boundaryChecker(xTry):
            self.updateBoundary()
            return xTry
        else:
            return self.x
    
    def getGradient(self, x):
        crx = []
        for j in range(len(self.wVec)):
            temp = 0.0
            for k in range(self.n):
                temp += self.tauMat[j][k] * self.x[k]
            crx.append(temp)
        grad = []
        for i in range(self.n):
            comp = 0.0
            for j in range(len(self.wVec)):
                comp += self.wVec[j] * self.tauMat[j][i] / crx[j]
            grad.append(comp)
        return grad
    
    def hillClimbing(self):
        self.getInitialX()
        self.enOpt = self.f(self.x)
        print(self.enOpt)
        T = 100
        while T > 0:
            print('Remains ' + str(T))
            T -= 1
            shrinkage = 0.01
            grad = self.normalize(self.getGradient(self.x))
            print('Grad = ' + str(grad))
            newx = []
            for i in range(self.n):
                newx.append(self.x[i] + grad[i] * shrinkage)
            newx = self.normalize(newx)
            while not self.boundaryChecker(newx):
                shrinkage *= 0.5
                for i in range(self.n):
                    if i == self.n - 1:
                        newx[i] = self.x[i]
                        continue 
                    newx[i] = self.x[i] + grad[i] * shrinkage
                newx = self.normalize(newx)
            
            enCur = self.f(newx)
            print('x  =   ' + str(newx))
            print('f(x) = ' + str(enCur))
            
            if enCur > self.enOpt:
                self.enOpt = enCur
                self.x = newx
            else:
                print('Unbelievable')
                return self.x
                self.x = newx
            
        print('Final x  =   ' + str(self.x))
        print('Final f(x) = ' + str(self.f(self.x)))
        return self.x
                   
    def calpsi(self):
        psi = []
        for g in range(len(self.length)):
            sa = 0.0
            sb = 0.0
            for i in range(len(self.length[g])):
                sa += self.x[self.mergeIndex[g][i][i]] / (self.length[g][i] - self.K + 1)
            for i in range(len(self.length[g])):
                for j in range(len(self.length[g])):
                    if i < j:
                        sb += self.x[self.mergeIndex[g][i][j]]/ (self.K - 1)
            for i in range(len(self.length[g])):
                temp = self.x[self.mergeIndex[g][i][i]] / (self.length[g][i] - self.K + 1)
                psi.append(temp / (sa - sb))   
        print(psi) 
        return psi
        
optimizer = OptimizeMLE()
optimizer.initializeSpecificParameters()
optimizer.hillClimbing()
optimizer.calpsi()

exit(0)
print('=================================================')
optimizer.x = optimizer.normalize([76, 96, 116, 136, 24, 0, 0, 24, 0, 24])
optimizer.hillClimbing()
optimizer.calpsi()

