import numpy as np
from numpy import linalg as LA
import sys

def print_mat(matrix):
    print("")
    for i in matrix:
        for j in i:
            print("%8.6f"%j,end="\t")
        print("")
    print("")
mat = [
[1.0000000 ,  0.2483624 ,  0.0000000 , -0.0000000 ,  0.0000000 ,  0.0630068  , 0.0630068 ,  0.0630068 ,  0.0630068],
[0.2483624 ,  1.0000000 ,  0.0000000 ,  0.0000000 ,  0.0000000 ,  0.4936348  , 0.4936348 ,  0.4936348 ,  0.4936348],
[0.0000000 ,  0.0000000 ,  1.0000000 ,  0.0000000 ,  0.0000000 ,  0.2707690  , 0.2707690 , -0.2707690 , -0.2707690],
[-0.0000000,   0.0000000,   0.0000000,   1.0000000,   0.0000000,  -0.2707690 ,  0.2707690,   0.2707690,  -0.270769],
[0.0000000 ,  0.0000000 ,  0.0000000 ,  0.0000000 ,  1.0000000 , -0.2707690  , 0.2707690 , -0.2707690 ,  0.2707690],
[0.0630068 ,  0.4936348 ,  0.2707690 , -0.2707690 , -0.2707690 ,  1.0000000  , 0.1714170 ,  0.1714170 ,  0.1714170],
[0.0630068 ,  0.4936348 ,  0.2707690 ,  0.2707690 ,  0.2707690 ,  0.1714170  , 1.0000000 ,  0.1714170 ,  0.1714170],
[0.0630068 ,  0.4936348 , -0.2707690 ,  0.2707690 , -0.2707690 ,  0.1714170  , 0.1714170 ,  1.0000000 ,  0.1714170],
[0.0630068 ,  0.4936348 , -0.2707690 , -0.2707690 ,  0.2707690 ,  0.1714170  , 0.1714170 ,  0.1714170 ,  1.0000000]]

matL = [
[1.0000000],
[0.2483624 ,  1.0000000],
[0.0000000 ,  0.0000000 ,  1.0000000 ],
[-0.0000000,   0.0000000,   0.0000000,   1.0000000],
[0.0000000 ,  0.0000000 ,  0.0000000 ,  0.0000000 ,  1.0000000],
[0.0630068 ,  0.4936348 ,  0.2707690 , -0.2707690 , -0.2707690 ,  1.0000000],
[0.0630068 ,  0.4936348 ,  0.2707690 ,  0.2707690 ,  0.2707690 ,  0.1714170  , 1.0000000],
[0.0630068 ,  0.4936348 , -0.2707690 ,  0.2707690 , -0.2707690 ,  0.1714170  , 0.1714170 ,  1.0000000],
[0.0630068 ,  0.4936348 , -0.2707690 , -0.2707690 ,  0.2707690 ,  0.1714170  , 0.1714170 ,  0.1714170 ,  1.0000000]]


mat2=[
[   1.0000000,   0.2367039,  -0.0000000,  -0.0000000,  -0.0000000,   0.0384056,   0.0384056],
[   0.2367039,   1.0000000,   0.0000000,  -0.0000000,   0.0000000,   0.3861388,   0.3861388],
[  -0.0000000,   0.0000000,   1.0000000,  -0.0000000,   0.0000000,   0.2684382,  -0.2684382],
[  -0.0000000,  -0.0000000,  -0.0000000,   1.0000000,  -0.0000000,   0.2097269,   0.2097269],
[  -0.0000000,   0.0000000,   0.0000000,  -0.0000000,   1.0000000,  -0.0000000,  -0.0000000],
[   0.0384056,   0.3861388,   0.2684382,   0.2097269,  -0.0000000,   1.0000000,   0.1817599],
[   0.0384056,   0.3861388,  -0.2684382,   0.2097269,  -0.0000000,   0.1817599,   1.0000000]]

w,v=LA.eigh(mat)
for i in range(len(w)):
    print(w[i])
    print(v[i])
v2=np.transpose(v)
vv=np.matmul(v2,v)
print_mat(vv)
w3,v3=LA.eig(vv)
print(w3)
print_mat(v3)
new_mat=np.zeros((len(v2),len(v2)))
for i in range(len(v2)):
    new_mat[i][i]=1/np.sqrt(w[i])
tempmat=np.matmul(np.matmul(v,new_mat),v2)
print_mat(tempmat)
