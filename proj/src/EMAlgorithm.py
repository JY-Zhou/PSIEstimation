import scipy as scp
import scipy.optimize as opt
import scipy.sparse as spa
import numpy as np
import json
from errno import EIDRM

class EMAlgorithm:
    def __init__(self, kmerHasher):
        self.NG = kmerHasher.NG
        self.NE = kmerHasher.NE
        self.K = kmerHasher.K
        self.NW = kmerHasher.NW
        self.readLength = kmerHasher.readLength
        
        self.NX = []
        for g in range(self.NG):
            self.NX.append(int(np.round(0.5 * (self.NE[g] * self.NE[g] + self.NE[g]))))
        self.NXSUM = [0]
        for g in range(self.NG):
            self.NXSUM.append(self.NX[g] + self.NXSUM[g])
        
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
        row = 0
        for kmer in kmerHasher.kmerTable:
            self.W[0, row] = kmerHasher.kmerTable[kmer][0]
            contribution = kmerHasher.kmerTable[kmer][1]
            for loc in contribution:
                sub = loc.split(',')
                g = int(sub[0])
                ei = int(sub[1])
                ej = int(sub[2])
                col = self.mergeIndex(g, ei, ej)                                
                self.Tau[row, col] = contribution[loc] / self.L[g][0, col - self.NXSUM[g]]
            row += 1
            
        for x in self.Tau:
            print(x)
            
        self.NA = []
        for g in range(self.NG):
            self.NA.append(int(np.round(0.5 * (self.NE[g] * self.NE[g] + 5 * self.NE[g] - 4))))
        self.A = []
        for g in range(self.NG):
            self.A.append(np.zeros((self.NA[g], self.NX[g])))
        for g in range(self.NG):
            row = 0
            
            #...Update A...
            
    def mergeIndex(self, g, ei, ej):
        e = 0
        if ei == ej:
            e = ei
        else:
            e = self.NE[g] + (2*self.NE[g] - ei - 3) * ei / 2 + ej - 1            
        return e + self.NXSUM[g]
    
    def splitIndex(self, col):
        #...Split Index...
        #return (g, ei, ej)
        return

    def initialVariables(self):
        self.Z = np.random.rand(1, self.NG)
        self.X = []
        for g in range(self.NG):
            self.X.append(np.random.rand(1, self.NX[g]))
        self.Mu = spa.lil_matrix((self.NW, self.NG))
        return
    
    def eStep(self):
        for s in range(self.NW):
            tot = 0.0
            for g in range(self.NG):
                self.Mu[s, g] = self.Z[0, g] * self.Tau[s,self.NXSUM[g]:self.NXSUM[g+1]].dot(self.X[g].T)[0, 0]
                tot += self.Mu[s, g]
            for g in range(self.NG):
                self.Mu[s, g] /= tot 
        return
    
    def mStep(self):
        tot = 0.0
        for g in range(self.NG):
            self.Z[0, g] = np.dot(self.W, self.Mu[:,g])
            tot += self.Z[0, g]
        for g in range(self.NG):
            self.Z[0, g] /= tot            
        for g in range(self.NG):
            self.optimizeQ(g)
        return
    
    def optimizeQ(self, g):
        #=======================================================================
        # print(np.shape(self.A[g]))
        # I = (np.dot(self.A[g], self.X[g].T) == 0.)
        # idxA1 = []
        # for r in range(self.NA[g]):
        #     if I[r, 0]:
        #         idxA1.append(r)
        # A1 = self.A[g][idxA1,:]
        # E = np.ones((1, self.NE[g]))
        # print(A1)
        # print(E)
        # input()
        #=======================================================================
        X0 = (g, self.X[g].T)
        cons = ({'type': 'ineq', 'fun': lambda X: np.dot(self.A[X[0]], X[1])}, 
                {'type': 'eq', 'fun': lambda X: np.sum(X[1]) - 1})
        #res = opt.minimi
        ze(self.objectFunction, x0 = X0, constraints = cons)
        #self.X[g] = res.x[1]
        return
    
    def objectFunction(self, X):
        res = 0.0
        for s in range(self.NW):
            res += self.W[s] * self.Mu[s][X[0]] * np.log(self.Tau[s][X[0]], X[1])
        return res
    
    def work(self, time):
        self.initialVariables()
        proc = 0
        while time > 0:
            if proc % 10 == 0:
                print(str(proc) + ' iteration processed...')
            proc += 1
            #print(self.Z) 
            self.eStep()
            self.mStep()
            time -= 1
        self.computePSI()
        return
    
    def computePSI(self):
        #=======================================================================
        # #Unit test for computing PSI
        # self.X = [np.matrix([[10,7,3,10,7,3,0,0,7,3]])]
        # self.L = [np.matrix([[1 ,1 ,1,1 ,1 ,1,1,1,1,1]])]
        # self.NG = 1
        # self.NX = [10]
        # self.NE = [4]
        #=======================================================================
               
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