# Problem 2
using LinearAlgebra
using Plots
using LaTeXStrings
include("readclassjson.jl");
print("\033c")

data = readclassjson("Midterm/quad_game.json");
b = data["b"]
x = data["xstar"]
beta = data["beta"]

A = zeros((4,4));
xx = [x[2] x[3] x[4] 0 0 0;x[1] 0 0 x[3] x[4] 0;0 x[1] 0 x[2] 0 x[4];0 0 x[1] 0 x[2] x[3]];
aa = xx\(b-2*beta.*x);
k = 1;
for i in 1:3
    for j in i+1:4
        A[i,j] = aa[k]
        aa
        A[j,i] = aa[k]
        k += 1
    end
end
print("A = ")
display(A)
println("\n(A+2βI)x* = $((A+2*beta*I)*x)")
println("b = $b")


