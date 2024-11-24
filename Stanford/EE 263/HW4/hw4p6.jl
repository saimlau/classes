# Problem 6
using LinearAlgebra
using Plots
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW4/piecepoly.json")
g = data["g"]
tau = data["tau"]
s = data["s"]
d = data["d"]

interval = [0,4];
k = length(tau)+1

# part b
function fit(ss,gs, d)
    m = length(ss)
    A = ones((m,d+1))
    for j in 1:d
        A[:,j+1] = ss.^j
    end
    return A\gs
end
# part d
function fit_c(ss,gs, d, g_offset, s_offset)
    m = length(ss)
    A = zeros((m,d))
    for j in 1:d
        A[:,j] = (ss.-s_offset).^j
    end
    return A\(gs.-g_offset)
end
scatter(s,g,label="data points")

coeff = zeros((k,d+1))
coeff_c = zeros((k,d+1))
s_tem1 = copy(s)
g_tem1 = copy(g)
cost = 0.0
cost_c = 0.0
for i in 1:k-1
    indx = findall(x->x<tau[i],s_tem1)
    s_tem = copy(s_tem1[indx])
    g_tem = copy(g_tem1[indx])
    deleteat!(s_tem1, indx)
    deleteat!(g_tem1, indx)
    coeff[i,:] = fit(s_tem,g_tem,d)
    if i==1
        coeff_c[i,:] = fit(s_tem,g_tem,d)
        cost_c += norm(g_tem-[evalpoly(si, coeff_c[i,:]) for si in s_tem])^2
    else
        coeff_c[i,1] = evalpoly(tau[i-1], coeff[i-1,:])
        coeff_c[i,2:end] = fit_c(s_tem,g_tem,d,coeff_c[i,1], tau[i-1])
        cost_c += norm(g_tem-[evalpoly(si-tau[i-1], coeff_c[i,:]) for si in s_tem])^2
    end
    cost += norm(g_tem-[evalpoly(si, coeff[i,:]) for si in s_tem])^2

    if i==1
        S = interval[1]:0.01:tau[i]
    else
        S = tau[i-1]:0.01:tau[i]
    end
    plot!(S,[evalpoly(si, coeff[i,:]) for si in S], label="p$i")
    if i==1
        plot!(S,[evalpoly(si, coeff_c[i,:]) for si in S], label="p$i continuous")
    else
        plot!(S,[evalpoly(si-tau[i-1], coeff_c[i,:]) for si in S], label="p$i continuous")
    end
end
coeff[k,:] = fit(s[s.>tau[k-1]],g[s.>tau[k-1]],d)
coeff_c[k,1] = evalpoly(tau[k-1], coeff[k-1,:])
coeff_c[k,2:end] = fit_c(s_tem1,g_tem1,d,coeff_c[k,1], tau[k-1])
S = tau[k-1]:0.01:interval[2]
plot!(S,[evalpoly(si, coeff[k,:]) for si in S], label="p$k")
plot!(S,[evalpoly(si-tau[k-1], coeff_c[k,:]) for si in S], label="p$k continuous")

savefig("HW4/piecePoly2.png")
println("Cost = $cost")
println("Coeff = ")
display(coeff)
println("Cost_continuous = $cost_c")


