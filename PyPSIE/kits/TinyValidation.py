import scipy as scp
import scipy.optimize as opt
import numpy as np

y = [1,1,0,1,0,0,1,0,1,1]

def likelihood(x):
    ret = 1
    for i in y:
        ret *= x[0] * (x[1]**i) * ((1-x[1])**(1-i)) + (1-x[0]) * (x[2]**i) * ((1-x[2]) ** (1-i))
    return -ret 

res = opt.minimize(fun = likelihood, 
                   x0 = [0.5, 0.5, 0.5],
                   bounds = [(0,1) for i in range(3)])
print(res)
res = opt.minimize(fun = likelihood, 
                   x0 = [0.4, 0.6, 0.7],
                   bounds = [(0,1) for i in range(3)])
print(res)