# Problem 4
using LinearAlgebra
using Plots
using LaTeXStrings
include("readclassjson.jl");

print("\033c")
data = readclassjson("Final/neural_network.json");
y = data["y"];
X = data["X"];
m = data["m"];
N = data["N"];
u = data["u"];
μ = data["mu"];
n = data["n"];

function AA(u)
    A = zeros(m,2*N)
    for i in 1:m
        for j in 1:N
            A[i,j*2-1:j*2] = -2(dot(X[i,:],u[j*2-1:j*2]))*X[i,:]
        end    
    end
    return A
end
function rr(u)
    r = []
    for i in 1:m
        push!(r, y[i]-sum([dot(X[i,:],u[j*2-1:j*2])^2 for j in 1:N]))        
    end
    return r
end
b(u) = -rr(u);
function duu(u)
    B = [AA(u);sqrt(μ)*I]
    d = vcat(b(u),zeros(N*n))
    du = inv(transpose(B)*B)*transpose(B)*d
    return du
end
K_max = 300;
ress = [norm(rr(u))];
final_k = K_max
for k in 1:K_max
    u += duu(u)
    push!(ress,norm(rr(u)))
    if k>2 && abs(ress[end]-ress[end-1])<1e-5
        final_k = k
        break
    end
end
plot(0:final_k,ress, label="||r(u⁽ᵏ⁾)||", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color = "purple")
xlabel!("k")
ylabel!("||residual||")
savefig("Final/p4c.png")