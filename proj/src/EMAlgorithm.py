import scipy as scp
import scipy.optimize as opt
import scipy.sparse as spa
import numpy as np
import json

class EMAlgorithm:
    def __init__(self, kmerHasher):
        self.EPS = 1e-30
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
            self.NX.append(int(self.NE[g] * (self.NE[g] + 1) / 2))
            
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
                self.L[g][0, col] = ed - st - self.readLength + 1
                col += 1
            for ei in range(self.NE[g]):
                ej = ei + 1
                while ej < self.NE[g]:
                    self.L[g][0, col] = 2 * self.readLength - 2 - self.readLength + 1
                    ej += 1
                    col += 1
        
        self.Tau = spa.lil_matrix((self.NW, self.NXSUM[self.NG]))
        self.W = np.zeros((1, self.NW))
        self.MuNonZero = []
        row = 0
        for kmer in kmerHasher.kmerTable:
            self.W[0, row] = kmerHasher.kmerTable[kmer][0]
            contribution = kmerHasher.kmerTable[kmer][1]
            
            gSet = {}
            for loc in contribution:
                sub = loc.split(',')
                g = int(sub[0])
                ei = int(sub[1])
                ej = int(sub[2])
                col = self.MergeIdx[(g, ei, ej)]
                self.Tau[row, col] = contribution[loc]
                gSet[g] = True
                
            for g in gSet:
                self.MuNonZero.append((row, g))
            row += 1
            
    def initialConstraints(self):
        self.NA = []
        for g in range(self.NG):
            if self.NE[g] > 2:
                self.NA.append(3 * self.NE[g] - 4)
            else:
                self.NA.append(self.NE[g])
            
        self.A = []
        for g in range(self.NG):
            self.A.append(np.zeros((self.NA[g], self.NX[g])))
        for g in range(self.NG):
            if self.NA[g] == 1:
                self.A[g][0, 0] = 1
                break
            
            row = 0
            
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
            
            NPsi = NLeftJunc + self.NE[g] - 2
            while row < NPsi:
                for i in range(self.NX[g]):
                    if i >= self.NE[g]:
                        self.A[g][row, i] = -1 
                    elif i != row - NLeftJunc + 1:
                        self.A[g][row, i] = 1
                row += 1                
                
            self.A[g] = self.A[g] / (np.ones((self.NA[g], 1)).dot(self.L[g]))
    
    def initialX(self, g):
        ret = np.random.rand(1, self.NX[g])
        tot = 0.0
        for e in range(self.NX[g]):
            tot += ret[0, e]
        for e in range(self.NX[g]):
            ret[0, e] /= tot
                
        while not (self.A[g].dot(ret.T) > -self.EPS).all():
            ret = np.random.rand(1, self.NX[g])
            for e in range(self.NE[g]):
                ret[0, e] *= 2
            tot = 0.0
            for e in range(self.NX[g]):
                tot += ret[0, e]
            for e in range(self.NX[g]):
                ret[0, e] /= tot
        return ret

    def initialVariables(self):
        self.Z = np.random.rand(1, self.NG)
        tot = 0.0
        for g in range(self.NG):
            tot += self.Z[0, g]
        for g in range(self.NG):
            self.Z[0, g] /= tot
        
        self.X = []
        for g in range(self.NG):
            self.X.append(self.initialX(g))

        self.Mu = spa.lil_matrix((self.NW, self.NG))
        return
    
    def eStep(self):
        tot = []
        for s in range(self.NW):
            tot.append(0)
        for loca in self.MuNonZero:
            s = loca[0]
            g = loca[1]
            self.Mu[s, g] = self.Z[0, g] * self.Tau[s, self.NXSUM[g]:self.NXSUM[g+1]].dot(self.X[g].T)[0, 0]
            tot[s] += self.Mu[s, g]
        for loca in self.MuNonZero:
            s = loca[0]
            g = loca[1]
            self.Mu[s, g] /= tot[s]
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
        glopt = float('inf')
        
        for t in range(15):
            xInit = self.initialX(g)
            print((self.A[0].dot(self.X[0].T) > 1e-15).all())
            print(np.ones((1, self.NX[0])).dot(self.X[0].T))
            res = opt.minimize(fun = self.QFunction,
                               x0 = xInit,
                               args = (g,),
                               tol = self.EPS, 
                               bounds = [(0, 1) for i in range(self.NX[g])],
                               #jac = self.QDerivate,
                               constraints = ({'type':'ineq', 'fun':lambda X: self.A[g].dot(X.T)},
                                              {'type':'eq', 'fun':lambda X: np.ones((1, self.NX[g])).dot(X.T) - 1}))
            print(res.fun)
            if res.fun[0, 0] < glopt:
                finres = res
                glopt = res.fun[0, 0]
        print(finres)
        for i in range(self.NX[g]):
            self.X[g] = np.matrix(finres.x)
        return
     
    def QFunction(self, X, g):
        X = np.matrix(X)
        temp = self.Tau[:,self.NXSUM[g]:self.NXSUM[g+1]].dot(X.T)
        if not (temp > self.EPS).all():
            return float('inf')
        return -self.Mu[:,g].multiply(np.log(temp)).T.dot(self.W.T)
    
    def QDerivate(self, X, g):
        X = np.matrix(X)
        den = self.Tau[:,self.NXSUM[g]:self.NXSUM[g+1]].dot(X.T)
        print(self.Mu[:,g].dot(np.ones((1, self.NX[g]))).shape)
        print(self.Tau[:,self.NXSUM[g]:self.NXSUM[g+1]].shape)
        print(den.shape)
        
        num = self.Mu[:,g].dot(np.ones((1, self.NX[g]))) * (self.Tau[:,self.NXSUM[g]:self.NXSUM[g+1]])
        print(num.shape)
        jac = self.W.dot(num) / (den.dot(np.ones((1, self.NX[g]))))
        return jac
     
    def likelihoodFunction(self, x):
        temp = np.zeros((self.NW, 1))
        x = np.matrix(x)
        for g in range(self.NG):
            temp += self.Z[0, g] * self.Tau[:,self.NXSUM[g]:self.NXSUM[g+1]].dot(x.T)
        temp = scp.log(temp)
        return -self.W.dot(temp)
     
    #===========================================================================
    # def optimizeLikelihood(self):
    #     res = opt.minimize(fun = self.likelihoodFunction,
    #                        x0 = self.X[0],
    #                        bounds = [(0, 1) for i in range(self.NX[0])],
    #                        constraints = ({'type':'ineq', 'fun': lambda x: self.A[0].dot(x.T)},
    #                                       {'type':'eq', 'fun': lambda x: np.ones((1, self.NX[0])).dot(x.T) - 1 }))
    #     return res
    #===========================================================================

    def work(self, time):
        self.initialVariables()
                
        #=======================================================================
        # res = self.optimizeLikelihood()        
        # print(res)
        # self.X[0] = np.matrix(res.x)
        #=======================================================================
               
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
            e = 0
            while e < self.NE[g]:
                sumEx += tempPsi[0, e]
                e += 1
            while e < self.NX[g]:
                sumJu += tempPsi[0, e]
                e += 1
            tempPsi /= (sumEx - sumJu)
            self.Psi.append(tempPsi[0,:self.NE[g]].tolist())
        print(self.Psi)
        psiFile = open('../output/PsiResult.json', 'w')
        json.dump(self.Psi, psiFile)

        return