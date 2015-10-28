import scipy as scp
import numpy as np
import json
from time import sleep

class EMAlgorithm:
    def __init__(self, NG, NE, NW, K):
        self.NG = NG
        self.NE = NE
        self.K = K
        self.NW = NW
        
        self.L = []
        for g in range(self.NG):
            self.L.append(np.random.rand(1, self.NE[g]))
        
        self.C = []
        for s in range(self.NW):
            self.C.append([])
            for g in range(self.NG):
                self.C[s].append(np.random.rand(1, self.NE[g]))
        
        self.W = np.random.rand(1, self.NW)
        self.Z = np.random.rand(1, self.NG)
        self.X = []
        for g in range(self.NG):
            self.X.append(np.random.rand(1, self.NE[g]))
            
        self.Mu = np.random.rand(self.NW, self.NG)
        self.A = []
        for g in range(self.NG):
            self.A = np.random.rand(self.numberOfConstraints(self.NE[g]), self.NE[g])

        self.initialByKmerTable()

    def numberOfConstraints(self, NEg):
        return np.round(0.5 *(NEg * NEg + 5 * NEg - 4))
    
    def initialByKmerTable(self):
        return
    
    def initialVariables(self):
        return        
    
    def eStep(self):
        for s in range(self.NW):
            tot = 0.0
            for g in range(self.NG):
                self.Mu[s, g] = self.Z[0, g] * np.dot(self.C[s][g]/self.L[g], self.X[g].T)[0, 0]
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
        
        return
          
    def work(self, time):
        self.initialByKmerTable()
        self.initialVariables()
        while time > 0:
            print(self.Z)
            self.eStep()
            self.mStep()
            time -= 1
        self.computePSI()
        return
    
    def computePSI(self):
        self.Psi = []
        for g in range(self.NG):
            self.Psi.append([])
            for e in range(self.NE[g]):
                self.Psi[g].append(0.0)
                    
        return
    
    def show(self):
        print(self.X)
        print(self.Z)
        
if __name__ == "__main__":
    test = EMAlgorithm(3, [2, 3, 5], 4, 25)
    test.work(100) 