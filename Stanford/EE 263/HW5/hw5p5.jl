# Problem 5
using LinearAlgebra
print("\033c")
G = [2 3;1 0;0 4;1 1;-1 2];
G_til = [-3 -1;-1 0;2 -3;-1 -3;1 2];
GG = vcat(G,G_til)

H = inv(transpose(GG)*GG)*transpose(GG);
display(H[:,1:5])
display(H[:,6:10])
