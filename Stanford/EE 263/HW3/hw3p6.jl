# Problem 6
using LinearAlgebra

print("\033c")
A = [1 2 1;1 -1 -2;-2 1 3;1 -1 -2;1 1 0]

m = size(A)[1]
n = size(A)[2]

Q, R = qr(A);

Q_1 = Matrix(Q);
r = rank(A)
Q_b = Q_1[:,1:r]

# B = rand(m-r,m)*10;
B = [1.0 0 0 1.0 0;0 1.0 0 0 1.0;0 0 1.0 1.0 0]

for i=1:m-r 
    for j=1:r 
        B[i,:] -= dot(B[i,:],Q_b[:,j])*Q_b[:,j]
    end
end

println(size(B))
B_disp = round.(B, sigdigits=5)
println("B =")
display("text/plain", B_disp)

# Test if it works
x = rand(n)*300;
y = A*x;
println("Test with random x in R^3:")
println(B*y)
println("Test when y inconsistent:")
y[3] += (rand(1)*20)[1]
y[5] -= (rand(1)*3)[1]
println(B*y)
