# Problem 3
using LinearAlgebra
using Plots
using LaTeXStrings
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW7/time_comp_data.json");
c = data["c"]
k = 1
n = length(c)

C = zeros(2*n-2*k-2,n);
m = 1;
for i in 2:n-k
    for j in 1:i-1
        if i-j>n || i-j<=0
            C[m,j] = 0
        else
            C[m,j] = c[i-j]
        end
    end
    m += 1
end
for i in n+2+k:2*n
    for j in 1:i-1
        if j>n
            continue
        else
            if i-j>n || i-j<=0
                C[m,j] = 0
            else
                C[m,j] = c[i-j]
            end
        end
    end
    m +=1
end
U, Σ, V = svd(C);
w = V[:,7];
fig = plot(1:n,w, label="w", layout=(1,2), size = (810, 360), framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color = "green", subplot=1)

h = zeros(2*n-1);
for i in 2:2*n
    for j in 1:i-1
        if !(j>n || i-j>n || j<=0 || i-j<=0)
            h[i-1] += w[j]*c[i-j]
        end
    end
end
Etot = norm(h)^2;
Edes = norm(h[n-k:n+k])^2;
DTE = Edes/Etot;
println("DTE = $DTE")

plot(fig, 2:2*n,h, label="h", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color = "red", subplot=2)
savefig("HW7/p3.png")

