# Problem 4
using LinearAlgebra
print("\033c")
A = [3 1;65 2;7 9;85 2];
B = [57 6;23 5;8 3;2 1];
x = [6; 7];
v = [3; 15];
w = [0.03;0.02;-0.1;0.05];
n = length(x);
p = length(v);

y = A*x+B*v+w;

x_N = A\y;
print("x_Nikola = ")
display(x_N)

xv_A = [A B]\y;
x_A = xv_A[1:n];
print("x_Almir = ")
display(x_A)
v_A = xv_A[n+1:n+p];

x_M = A\(y-B*v_A);
print("x_Miki = ")
display(x_M)