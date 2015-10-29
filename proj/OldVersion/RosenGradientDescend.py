import numpy

class RosenGradientDescend:
    EPS = 1e-12
    #Parameters of standard linear constraints
    A = numpy.zeros((0,0), dtype = numpy.float64)
    b = numpy.zeros((0,0), dtype = numpy.float64)
    E = numpy.zeros((0,0), dtype = numpy.float64)
    e = numpy.zeros((0,0), dtype = numpy.float64)
    N = 0
    
    def f(self, x):
        pass

    def gradient(self, x):
        pass
    
    def initialize(self):
        pass
    
    def findLamda(self, xk, dk, lamdaMx):
        pass
    
    def compute(self, x0):
        self.initialize()
        xk = x0
        while True:
            grad = self.gradient(xk)       
            #print('grad')
            #print(grad)     
            idxA1 = []
            idxA2 = []
            for i in range(self.A.shape[0]):
                if numpy.absolute(self.A[i,:].dot(xk) - self.b[i,:]) < self.EPS:
                    idxA1.append(i)
                else:
                    idxA2.append(i)
            A1 = self.A[idxA1,:]
            A2 = self.A[idxA2,:]
            b2 = self.b[idxA2,:]
            
            #print('A1,A2,b2')
            #print(A1)
            #print(A2)
            #print(b2)
            
            while True:
                if self.E.shape[0] == 0:
                    M = A1
                else:
                    M = numpy.concatenate((A1, self.E))
                    
                #print('M')
                #print(M)
                
                P = numpy.eye(self.N, dtype = numpy.float64)
                if M.shape[0] != 0:
                    #print('MMT')
                    mid = M.dot(M.T)
                    if numpy.linalg.matrix_rank(mid) < M.shape[0]:
                        print('End-3: M*M\' is Singular')
                        return xk
                    temp = numpy.linalg.inv(M.dot(M.T)).dot(M)
                    P = P - M.T.dot(temp)
                
                #print('P')
                #print(P)
                
                dk = -P.dot(grad)
                #print('dk')
                #print(dk)
                if (numpy.absolute(dk) < self.EPS).all():
                    if M.shape[0] == 0:
                        print('End-1: Inner optimal founded')
                        return xk
                    else:
                        W = temp.dot(grad)
                        u = W[:A1.shape[0],:]
                        
                        #print('u')
                        #print(u)
                        
                        if (u > -self.EPS).all():
                            print('End-2: K-T point founded')
                            return xk
                        else:
                            idxNew = []
                            flag = True
                            for i in range(u.shape[0]):
                                if flag and u[i, 0] < -self.EPS:
                                    flag = False
                                    continue
                                else:
                                    idxNew.append(i)
                            A1 = A1[idxNew,:]
                else:
                    break
            bCap = b2 - A2.dot(xk)
            dCap = A2.dot(dk)
            
            #print('bCap')
            #print(bCap)
            #print('dCap')
            #print(dCap)
            
            if (dCap > -self.EPS).all():
                lamdaMX = float('inf')
            else:
                lamdaMX = (bCap[(dCap < -self.EPS)] / dCap[(dCap < -self.EPS)]).min()
            lamda = self.findLamda(xk, dk, lamdaMX)
            print('lamdaMax, lamda')
            print(str(lamdaMX) + ', ' + str(lamda))
            xk = xk + lamda*dk
            #print('xk')
            #print(xk)
            #print('====================')
        
class OptimizeExample(RosenGradientDescend):    
    def f(self, x):
        return 2*x[0,0]*x[0,0] + 2*x[1,0]*x[1,0] - 2*x[0,0]*x[1,0] - 4*x[0,0] - 6*x[1,0]
    
    def gradient(self, x):
        grad = numpy.zeros((self.N, 1), dtype = numpy.float64)
        grad[0, 0] = 4*x[0, 0] - 2*x[1, 0] - 4
        grad[1, 0] = 4*x[1, 0] - 2*x[0, 0] - 6
        return grad

    def initialize(self):
        self.A = numpy.matrix([[-1, -1], 
                               [-1, -5],
                               [1, 0],
                               [0, 1]], dtype = numpy.float64)
        self.b = numpy.matrix([[-2],
                               [-5],
                               [0],
                               [0]], dtype = numpy.float64)
        self.N = 2
    
    def findLamda(self, xk, dk, lamdaMx):
        opti = - (self.gradient(xk).T.dot(dk))
        opti = opti / ((self.gradient(dk) + numpy.matrix([[4, 6]], dtype = numpy.float64).T).T.dot(dk))
        if opti[0,0] - lamdaMx > self.EPS:
            return lamdaMx
        if opti[0,0] < -self.EPS:
            return 0
        return opti[0,0]
                 
class OptimizeFirst(RosenGradientDescend):
    def f(self, x):
        return x[0,0]*x[0,0] + x[1,0]*x[1,0] + 2*x[1,0] + 5
    
    def gradient(self, x):
        grad = numpy.zeros((self.N, 1), dtype = numpy.float64)
        grad[0, 0] = 2*x[0, 0]
        grad[1, 0] = 2*x[1, 0] + 2
        return grad
    
    def initialize(self):
        self.A = numpy.matrix([[1, -2],
                              [1, 0],
                              [0, 1]], dtype = numpy.float64)
        self.b = numpy.zeros((3,1), dtype = numpy.float64)
        self.A = numpy.matrix([[1, -2],
                               [0, 1]], dtype = numpy.float64)
        self.b = numpy.zeros((2,1), dtype = numpy.float64)
        self.N = 2
    
    def findLamda(self, xk, dk, lamdaMx):
        opti = - (self.gradient(xk).T.dot(dk)) / (2*(dk[0,0]*dk[0,0] + dk[1,0]*dk[1,0]))
        if opti[0,0] - lamdaMx > self.EPS:
            return lamdaMx
        if opti[0,0] < -self.EPS:
            return 0
        return opti[0,0]

class OptimizeSecond(RosenGradientDescend):
    def gradient(self, x):
        grad = numpy.zeros((self.N, 1), dtype = numpy.float64)
        grad[0, 0] = 2*x[0, 0] + x[1,0] - 6
        grad[1, 0] = x[0, 0] + 4*x[1, 0] - 2
        grad[2, 0] = -12
        return grad
        
    def initialize(self):
        self.A = numpy.matrix([[1, -2, 0],
                               [1, 0, 0],
                               [0, 1, 0],
                               [0, 0, 1]], dtype = numpy.float64)
        self.b = numpy.matrix([[-3],
                               [0],
                               [0],
                               [0]], dtype = numpy.float64)
        self.E = numpy.matrix([[1, 1, 1]], dtype = numpy.float64)
        self.e = numpy.matrix([[2]], dtype = numpy.float64)
        self.N = 3
    
    def findLamda(self, xk, dk, lamdaMx):
        opti = - (self.gradient(xk).T.dot(dk))
        opti = opti / (2*dk[0,0]*dk[0,0] + 2*dk[0,0]*dk[1,0] + 4*dk[1,0]*dk[1,0])
        if opti[0,0] - lamdaMx > self.EPS:
            return lamdaMx
        if opti[0,0] < -self.EPS:
            return 0
        return opti[0,0]
    
optimizer = OptimizeExample()
res = optimizer.compute(numpy.matrix([[0, 0]], dtype = numpy.float64).T)
print('========RESULT========')
print(res)
optimizer = OptimizeFirst()
res = optimizer.compute(numpy.matrix([[2, 0]], dtype = numpy.float64).T)
print('========RESULT========')
print(res)
optimizer = OptimizeSecond()
res = optimizer.compute(numpy.matrix([[1, 0, 1]], dtype = numpy.float64).T)
print('========RESULT========')
print(res)
#===============================================================================
# a = numpy.matrix([[1,2],[3,4]])
# print(a[1,:])
# print(a[[],:].shape)
# print((a>1).all())
#===============================================================================
