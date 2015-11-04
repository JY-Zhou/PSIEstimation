import scipy as scp
import scipy.optimize as opt
import scipy.sparse as spa
import numpy as np
import json

class EMAlgorithm:
    def __init__(self, kmerHasher):
        self.NG = kmerHasher.NG
        self.NE = kmerHasher.NE
        self.K = kmerHasher.K
        self.NW = kmerHasher.NW
        self.readLength = kmerHasher.readLength
        
        self.initialIndice()
        self.initialCoefficients(kmerHasher)
        self.initialConstraints()
        
    def initialIndice(self):   
        self.NX = []
        for g in range(self.NG):
            self.NX.append(int(np.round(0.5 * (self.NE[g] * self.NE[g] + self.NE[g]))))
        self.NXSUM = [0]
        for g in range(self.NG):
            self.NXSUM.append(self.NX[g] + self.NXSUM[g])
            
        self.MergeIdx = {}
        self.SplitIdx = []
        idx = 0
        for g in range(self.NG):
            for e in range(self.NE[g]):
                self.SplitIdx.append((g,e,e))
                self.MergeIdx[(g,e,e)] = idx
                idx += 1
            for ei in range(self.NE[g]):
                for ej in range(ei + 1, self.NE[g]):
                    self.SplitIdx.append((g,ei,ej))
                    self.MergeIdx[(g,ei,ej)] = idx
                    idx += 1                 

    def initialCoefficients(self, kmerHasher):
        self.L = []
        for g in range(self.NG):
            self.L.append(np.zeros((1, self.NX[g])))
        for g in range(self.NG):
            col = 0
            for e in range(self.NE[g]):
                st = kmerHasher.geneBoundary[g][e][0]
                ed = kmerHasher.geneBoundary[g][e][1] + 1
                self.L[g][0, col] = ed - st - self.K + 1
                col += 1
            for ei in range(self.NE[g]):
                ej = ei + 1
                while ej < self.NE[g]:
                    self.L[g][0, col] = 2 * self.readLength - 2 - self.K + 1
                    ej += 1
                    col += 1
        
        self.Tau = spa.lil_matrix((self.NW, self.NXSUM[self.NG]))
        self.W = np.zeros((1, self.NW))
        self.MuNonZero = []
        row = 0
        for kmer in kmerHasher.kmerTable:
            self.W[0, row] = kmerHasher.kmerTable[kmer][0]
            contribution = kmerHasher.kmerTable[kmer][1]
            for loc in contribution:
                sub = loc.split(',')
                g = int(sub[0])
                ei = int(sub[1])
                ej = int(sub[2])
                col = self.MergeIdx[(g, ei, ej)]
                self.Tau[row, col] = contribution[loc] / self.L[g][0, col - self.NXSUM[g]]
                self.MuNonZero.append((row, g))
            row += 1
            
    def initialConstraints(self):
        self.NA = []
        for g in range(self.NG):
            self.NA.append(int(np.round(0.5 * (self.NE[g] * self.NE[g] + 5 * self.NE[g] - 4))))
            
        self.A = []
        for g in range(self.NG):
            self.A.append(np.zeros((self.NA[g], self.NX[g])))
        for g in range(self.NG):
            row = 0            
            #
            NRightJunc = self.NE[g] - 1
            l = self.NE[g]
            r = l + self.NE[g] - 1
            while row < NRightJunc:
                self.A[g][row, row] = 1
                for i in range(l, r):
                    self.A[g][row, i] = -1
                l = r
                r = r + self.NE[g] - (row + 2)
                row += 1
            #
            NLeftJunc = NRightJunc + self.NE[g] - 1
            while row < NLeftJunc:
                self.A[g][row, row - NRightJunc + 1] = 1
                l = self.NE[g]
                r = row - NRightJunc
                i = 1
                while r >= 0:
                    self.A[g][row, l + r] = -1
                    l += self.NE[g] - i
                    i += 1
                    r -= 1     
                row += 1                
            #
            NPsi = NLeftJunc + self.NE[g]
            while row < NPsi:
                for i in range(self.NX[g]):
                    if i >= self.NE[g]:
                        self.A[g][row, i] = -1 
                    elif i != row - NLeftJunc:
                        self.A[g][row, i] = 1
                row += 1                
            #
            NBeta = NPsi + self.NE[g] * (self.NE[g] - 1) / 2
            while row < NBeta:
                self.A[g][row, row - NPsi + self.NE[g]] = 1
                row += 1
        #print(self.A[0])

    def initialVariables(self):
        self.Z = np.random.rand(1, self.NG)
        tot = 0.0
        for g in range(self.NG):
            tot += self.Z[0, g]
        for g in range(self.NG):
            self.Z[0, g] /= tot
        
        self.X = []
        for g in range(self.NG):
            self.X.append(np.random.rand(self.NX[g], 1))
            tot = 0.0
            for e in range(self.NX[g]):
                tot += self.X[g][e, 0]
            for e in range(self.NX[g]):
                self.X[g][e, 0] /= tot
                            
            
        self.Mu = spa.lil_matrix((self.NW, self.NG))
        return
    
    def eStep(self):
        tot = []
        for g in range(self.NG):
            tot.append(0)
        for loca in self.MuNonZero:
            s = loca[0]
            g = loca[1]
            self.Mu[s, g] = self.Z[0, g] * self.Tau[s, self.NXSUM[g]:self.NXSUM[g+1]].dot(self.X[g])[0, 0]
            tot[g] += self.Mu[s, g]
        for loca in self.MuNonZero:
            s = loca[0]
            g = loca[1]
            self.Mu[s, g] /= tot[g]
        return
    
    def mStep(self):
        tot = 0.0
        for g in range(self.NG):
            self.Z[0, g] = self.Mu[:,g].T.dot(self.W.T)
            tot += self.Z[0, g]
        for g in range(self.NG):
            self.Z[0, g] /= tot

        for g in range(self.NG):
            self.optimizeQ(g)
        return
    
    def optimizeQ(self, g):
        print(self.Mu[:,g].multiply(scp.log(self.Tau[:,self.NXSUM[g]:self.NXSUM[g+1]].dot(self.X[g]))).T.dot(self.W.T)[0,0])
        res = opt.minimize(fun = lambda X: -self.Mu[:,g].multiply(scp.log(self.Tau[:,self.NXSUM[g]:self.NXSUM[g+1]].dot(X))).T.dot(self.W.T)[0,0],
                           x0 = self.X[g],
                           bounds = [(0, 1) for i in range(self.NX[g])],
                           constraints = ({'type':'ineq', 'fun':lambda X: self.A[g].dot(X)},
                                          {'type':'eq', 'fun':lambda X: np.ones((1, self.NX[g])).dot(X) - 1}))

        print(res)
        input()
        for i in range(self.NX[g]):
            self.X[g][i, 0] = res.x[i]
        return
    
    def logarithm(self, X):
        return         
    
    def work(self, time):
        self.initialVariables()
        proc = 0
        while proc < time:
            if proc % 1 == 0:
                print(str(proc) + ' iteration processed...')
            proc += 1
            self.eStep()
            self.mStep()
        self.computePSI()
        return
    
    def computePSI(self):
        self.Psi = []
        for g in range(self.NG):
            tempPsi = self.X[g] / self.L[g]  
            sumEx = 0.0
            sumJu = 0.0
            for e in range(self.NX[g]):
                if e < self.NE[g]:
                    sumEx += tempPsi[0, e]
                else:
                    sumJu += tempPsi[0, e]
                e += 1
            tempPsi = tempPsi / (sumEx - sumJu)
            self.Psi.append(tempPsi[0,:self.NE[g]].tolist()) 
        psiFile = open('../output/PsiResult.json', 'w')
        json.dump(self.Psi, psiFile)
        return