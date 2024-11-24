# Problem 6
using LinearAlgebra
using LaTeXStrings
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW7/regl_data.json");
n = data["n"]
c = data["c"]
m = data["m"]
w = data["w"]
x = data["x"]
m1 = data["m1"]
h = data["h"]
d = data["d"]

# part a
A = zeros(n,n);
for i in 1:n
    for k in -h:h
        if i+k<1 || i+k>n
            continue
        end
        A[i,i+k] = c[k+h+1]
    end
end
U,Σ,V = svd(A);
r = length(Σ);
plot(1:r, Σ, label="σₖ", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color = "brown")
savefig("HW7/p6a.png")

# part b
plot(1:200, V[:,1], label="v1", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottomright, color = "red")
colors = ["red","green","blue","purple","brown"];
for j in 2:6
    plot!(1:200, V[:,j], label="v$j", color = colors[j-1], linewidth=2)
end
savefig("HW7/p6b.png")

# part c
ymeas = A*x+w;
xls = V*inv(Diagonal(Σ))*transpose(U)*ymeas;
plot(1:n, xls, label="xₗₛ", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color = "brown")
savefig("HW7/p6c.png")

# part d 
plot()
for r in [5,10,15,30,50]
    xls = V[:,1:r]*inv(Diagonal(Σ[1:r]))*transpose(U[:,1:r])*ymeas;
    plot!(1:n, xls, label="r=$r", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
            linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright)
end
savefig("HW7/p6d.png")

# part e
err = []
for r in 1:35
    xls = V[:,1:r]*inv(Diagonal(Σ[1:r]))*transpose(U[:,1:r])*ymeas;
    push!(err, norm(x-xls))
end
plot(1:35, err, label="error", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright)
savefig("HW7/p6e.png")

# part f
r = 21;
xls = V[:,1:r]*inv(Diagonal(Σ[1:r]))*transpose(U[:,1:r])*ymeas;
plot(1:n, xls, label="r=$r", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright)
savefig("HW7/p6f.png")

# part g
err = []
mus = 0.01:0.01:0.5;
for mu=mus
    B = [A;sqrt(mu)*I(n)]
    yy = vcat(ymeas,zeros(n))
    UU,SS,VV = svd(B)
    xreg = VV*inv(Diagonal(SS))*transpose(UU)*yy;
    push!(err, norm(x-xreg))
end
plot(mus, err, label="error", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright)
savefig("HW7/choosingMu.png")
mu = 0.05;
B = [A;sqrt(mu)*I(n)];
yy = vcat(ymeas,zeros(n));
UU,SS,VV = svd(B);
xreg = VV*inv(Diagonal(SS))*transpose(UU)*yy;
plot(1:n, xreg, label="μ=$mu", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright)
savefig("HW7/p6g.png")

