# Problem 5 c
using LinearAlgebra
using LaTeXStrings
import PlotlyJS
include("readclassjson.jl");
print("\033c")
data = readclassjson("Midterm/brain_data_2.json")
n = data["n"]
U = data["U"]
k = data["k"]
T = data["T"]
X = Dict()
X[1] = data["X_1"]
X[2] = data["X_2"]
X[3] = data["X_3"]

function xxx(i,t,n,X)
    xx = zeros((n,n^2))
    for j in 1:n
        xx[j,(j-1)*n+1:j*n] = X[i][:,t-1];
    end
    return xx
end

XX = zeros((0,n^2));
b = zeros((0,1));
for i in 1:k
    for t in 2:T
        XX = vcat(XX,xxx(i,t,n,X))
    end
    b = vcat(b,vec(X[i][:,2:T]))
end

vW = inv(transpose(XX)*XX)*transpose(XX)*b;
W = transpose(reshape(vW,(n,n)));
println("Estimated W = ")
display(W)
PlotlyJS.savefig(PlotlyJS.plot(PlotlyJS.heatmap(z=W)),"Midterm/W2_est.png")