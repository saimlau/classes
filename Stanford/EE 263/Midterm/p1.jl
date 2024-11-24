using LinearAlgebra

A = [1 0;0 1;0 0]
x = [3,2,1]

transpose(A)*x

A = [1 0;0 1]
transpose(A)*A
B = [1 0;0 1;0 0]
transpose(B)*B