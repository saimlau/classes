import numpy as np
import array_to_latex as a2l

S = np.array([[1,1,1],[1,1,1],[1,1,1]])/9
D = np.array([[-1,0,1],[-2,0,2],[-1,0,1]])
A = np.array([[0,9,0,9,0],[0,9,0,9,0],[0,9,0,9,0],[0,9,0,9,0],[0,9,0,9,0]])
N = 1
M = 1

A_new = np.zeros((9,9))
A_new[2:7,2:7] = A
resl1 = np.zeros((7,7))
for i in range(7):
    for j in range(7):
        for k in range(3):
            for r in range(3):
                resl1[i,j] += S[k,r]*A_new[i+2-k,j+2-r]
print(resl1)
resl2 = np.zeros((5,5))
for i in range(5):
    for j in range(5):
        for k in range(3):
            for r in range(3):
                resl2[i,j] += D[k,r]*resl1[i+2-k,j+2-r]
print(resl2)
a2l.to_ltx(resl2, frmt = '{:6.0f}', arraytype = 'bmatrix', print_out=True)