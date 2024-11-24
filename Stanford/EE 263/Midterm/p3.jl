# Problem 3
using LinearAlgebra
using Plots
using LaTeXStrings
using Colors
include("readclassjson.jl");
print("\033c")

data = readclassjson("Midterm/ImageReconstructionData.json");
Y = data["Y"]
noi_Y = data["noised_Y"]

# Part c
plot(Gray.(Y));
title!("Given Y");
savefig("Midterm/orig_Y.png");
W = [1 -2 1;1 -2 1;1 -2 1];
mn1 = size(W)[1];
m = size(Y)[1];
n = m-mn1+1;
vY = vec(transpose(Y));
Ω = zeros((m,n,mn1));
D = zeros((m^2,n^2));
OO = zeros((mn1*m,n));
for i in 1:mn1
    for j in 1:n
        Ω[j:j+mn1-1,j,i] = W[i,:];
    end
    OO[(i-1)*m+1:(i-1)*m+m,1:n] = Ω[:,:,i]HW5
end
for i in 1:n
    D[(i-1)*m+1:(i-1)*m+mn1*m,(i-1)*n+1:(i-1)*n+n] = OO
end
vX = D\vY;
X = transpose(reshape(vX,(n,n)));
plot(Gray.(X));
title!("Recovered X");
savefig("Midterm/recovered_X.png");

# Part d
plot(Gray.(noi_Y));
title!("Given noised_Y");
savefig("Midterm/orig_noised_Y.png");
W = [2 -2 2;-2 0 -2;2 -2 2];
mn1 = size(W)[1];
m = size(noi_Y)[1];
n = m-mn1+1;
vY = vec(transpose(noi_Y));
Ω = zeros((m,n,mn1));
D = zeros((m^2,n^2));
OO = zeros((mn1*m,n));
d = ones((2*n+32)).*0.5;
for i in 1:mn1
    for j in 1:n
        Ω[j:j+mn1-1,j,i] = W[i,:];
    end
    OO[(i-1)*m+1:(i-1)*m+m,1:n] = Ω[:,:,i]
end
for i in 1:n
    D[(i-1)*m+1:(i-1)*m+mn1*m,(i-1)*n+1:(i-1)*n+n] = OO
end
C = zeros((2*n+32,n^2));
C[1:n,6*n+1:7*n] = Matrix(I,(n,n));
C[n+1:2*n,23*n+1:24*n] = Matrix(I,(n,n));
for i in 1:16
    C[2*n+i,7*n+(i-1)*n+1] = 1
    C[2*n+16+i,7*n+i*n] = 1
end
DT = transpose(D);
DTD = DT*D;
DTD_inv = inv(DTD);
CT = transpose(C);
vX = DTD_inv*(DT*vY-CT*inv(C*DTD_inv*CT)*(C*DTD_inv*DT*vY-d));
X = transpose(reshape(vX,(n,n)));
plot(Gray.(X));
title!("Recovered noised X");
savefig("Midterm/recovered_noised_X.png");
