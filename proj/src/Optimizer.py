import numpy

from src.RosenGradientDescend import RosenGradientDescend

class Optimizer(RosenGradientDescend):
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