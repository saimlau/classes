# Problem 5
using LinearAlgebra
using LaTeXStrings
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW7/opt_bin_data.json");
A = data["A"]
Vmax = data["Vmax"]
N = data["N"]

U,Σ,V = svd(A);
c = 2*Vmax/norm(Σ[2]*U[:,2]-Σ[1]*U[:,1]);
s1 = c*V[:,1];
s2 = c*V[:,2];
println("s₁=$(round.(s1,sigdigits=6))")
println("s₂=$(round.(s2,sigdigits=6))")