import numpy as np
import scipy.linalg as linalg

e = 1
A = np.zeros((3*e-2, e*(e+1)/2))
B = np.zeros((3*e-4, e*(e+1)/2))
row = 0            
NRightJunc = e - 1
l = e
r = l + e - 1
while row < NRightJunc:
    A[row, row] = 1
    B[row, row] = 1
    for i in range(l, r):
        A[row, i] = -1
        B[row, i] = -1
    l = r
    r = r + e - (row + 2)
    row += 1
    
NLeftJunc = NRightJunc + e - 1
while row < NLeftJunc:
    A[row, row - NRightJunc + 1] = 1
    B[row, row - NRightJunc + 1] = 1
    l = e
    r = row - NRightJunc
    i = 1
    while r >= 0:
        A[row, l + r] = -1
        B[row, l + r] = -1
        l += e - i
        i += 1
        r -= 1     
    row += 1           
    
NPsi = NLeftJunc + e
while row < NPsi:
    for i in range(int(e*(e+1)/2)):
        if i >= e:
            A[row, i] = -1 
        elif i != row - NLeftJunc:
            A[row, i] = 1
    row += 1

row = NLeftJunc
NPsi = NLeftJunc + e - 2
while row < NPsi:
    for i in range(int(e*(e+1)/2)):
        if i >= e:
            B[row, i] = -1 
        elif i != row - NLeftJunc + 1:
            B[row, i] = 1
    row += 1
        
print('---A---')     
print(A)
print(np.linalg.matrix_rank(A))
print('---B---')
print(B)
print(np.linalg.matrix_rank(B))

print('---Comp---')
print(linalg.lu(A)[-1])
print('')
print(linalg.lu(B)[-1])
print('')
l = []
for i in range(linalg.lu(A)[-1].shape[0]):
    if not (linalg.lu(A)[-1][i, :] == 0.).all():
        l.append(i)
d = linalg.lu(A)[-1][l, :] - linalg.lu(B)[-1]
print((d == 0.).all())
