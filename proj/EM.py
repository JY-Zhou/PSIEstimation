import scipy as scp
import numpy as np
import json

class EMAlgorithm:
    def __init__(self, NG, NE, NW, K):
        self.NG = NG
        self.NE = NE
        self.K = K
        self.NW = NW
        
        self.L = []
        for g in range(self.NG):
            self.L.append(np.zeros((1, self.NE[g])))
        
        self.C = []
        for s in range(self.NW):
            self.C.append([])
            for g in range(self.NG):
                self.C[s].append(np.zeros((1, self.NE[g])))
        
        self.W = np.zeros((1, self.NW))
        self.Z = np.zeros((1, self.NG))
        self.X = []
        for g in range(self.NG):
            self.X.append(np.zeros((1, self.NE[g])))
            
        self.Mu = np.zeros((self.NW, self.NG))
        self.A = []
        for g in range(self.NG):
            self.A = np.zeros((self.numberOfConstraints(self.NE[g]), self.NE[g]))

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
                self.Mu[s, g] = self.Z[g, 1] * np.dot(self.X[g], self.C[s][g]/self.L[g])
                tot += self.Mu[s, g]
            for g in range(self.NG):
                self.Mu[s, g] /= tot
        return
    
    def mStep(self):
        tot = 0.0
        for g in range(self.NG):
            self.Z[g, 1] = np.dot(self.W, self.Mu[:,g])
            tot += self.Z[g, 1]
        for g in range(self.NG):
            self.Z[g, 1] /= tot
            
        for g in range(self.NG):
            self.OptimizeQg()
        
        return
    
    def optimizeQg(self):
        return
          
    def work(self, time):
        self.initialByKmerTable()
        self.initialVariables()
        while time > 0:
            self.eStep()
            self.mStep()
            time -= 1
        return
        
if __name__ == "__main__":
    test = EMAlgorithm(3, [2, 3, 5], 57, 25)
    print(test.L)
    test.work(1)
    print(test.L)