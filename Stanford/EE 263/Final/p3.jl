# Problem 3
using LinearAlgebra
using Plots
print("\033c")

A = I + diagm(1=>[1.,0.,1.]);
B = [.5 0;1 0;0 .5;0 1];
Bg = [0,0,-1/2,-1];

# Part a
function CC(T)
    CT = zeros(4,0)
    for k in 0:T-1
        CT = hcat(CT,A^k*B)
    end
    return CT
end
function DD(T)
    DT = zeros(4)
    for τ in 0:T-1
        DT += A^τ*Bg
    end
    return DT
end
CCh(T) = CC(T)[[1,3],:];
DDh(T) = DD(T)[[1,3],:];
p=[2 2 0 -2 -2 0 0 0 0 0 0; 1/2 -2 0 1/2 -2 0 0 0 0 0 0];
T=[10 20 30 40 50 60 61 62 63 64 65];
N=length(T);
Â = zeros(N*2,2*65);
tem = zeros(N*2);
for i in 1:N
    Â[i*2-1:i*2,:] = [zeros(2,(65-T[i])*2) CCh(T[i])]
    tem[i*2-1:i*2] = DDh(T[i])
end
ŷ = vec(p)-tem;
useq = transpose(Â)*inv(Â*transpose(Â))*ŷ;
ud = reshape(useq,2,:);
ud = reverse(ud,dims=2);
xx = zeros(4,T[end]+1);
for t in 1:T[end]
    xx[:,t+1] = A*xx[:,t]+B*ud[:,t]+Bg
end
plot(xx[1,:], xx[3,:], label="trajectory", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottom, color = "green")
scatter!(p[1,:],p[2,:], markershape=:star5, label="way-points", color="purple")
xlabel!("x₁")
ylabel!("x₃")
savefig("Final/p3a1")
plot(0:T[end]-1, ud[1,:], label="u₁", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:right, color = "red")
plot!(0:T[end]-1,ud[2,:], label="u₂", linewidth=2, color="blue")
xlabel!("time")
ylabel!("Thrust")
savefig("Final/p3a2")

# Part b
function uuls(μ)
    Ã = [Â;sqrt(μ)*I]
    ỹ = vcat(ŷ,zeros(2*T[end]))
    useqq = inv(transpose(Ã)*Ã)*transpose(Ã)*ỹ
    udd = reshape(useqq,2,:)
    udd = reverse(udd,dims=2)
    J1 = norm(Â*useqq-ŷ)^2
    J2 = norm(useqq)^2
    return udd, J1, J2
end
function simulate(udd)
    xxs = zeros(4,T[end]+1);
    for t in 1:T[end]
        xxs[:,t+1] = A*xxs[:,t]+B*udd[:,t]+Bg
    end
    return xxs
end
μs = 10.0.^(-5:0.1:20);
J1s = [];
J2s = [];
for μ in μs
    ud, J1, J2 = uuls(μ)
    push!(J1s, J1)
    push!(J2s, J2)
end
plot(J2s, J1s, label="μ∈[10⁻⁵,10²⁰]", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color = "brown")
xlabel!("J₂")
ylabel!("J₁")
savefig("Final/p3b.png")

# Part c
plot()
for μ in [0.01, 0.1, 1, 10, 20]
    ud, J1, J2 = uuls(μ)
    xx = simulate(ud)
    plot!(xx[1,:], xx[3,:], label="μ=$μ", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottom)
end
scatter!(p[1,:],p[2,:], markershape=:star5, label="way-points", color="purple")
xlabel!("x₁")
ylabel!("x₃")
savefig("Final/p3c")

# Part f
B̂= zeros(65*2,2*65);
d = zeros(65*2);
for i in 1:65
    B̂[i*2-1:i*2,:] = [zeros(2,(65-i)*2) CCh(i)]
    d[i*2-1:i*2] = DDh(i)
end
function uucmols(γ)
    F = [B̂;sqrt(γ)*I]
    G = vcat(d, zeros(2*65))
    Σ = transpose(F)*F
    iS = inv(Σ)
    uopt = iS*transpose(Â)*inv(Â*iS*transpose(Â))*(ŷ+Â*iS*transpose(F)*G)-iS*transpose(F)*G
    J2 = norm(uopt)^2
    J3 = norm(B̂*uopt+d)^2
    udd = reshape(uopt,2,:)
    udd = reverse(udd,dims=2)
    return udd, J2, J3
end
γs = [1,10,100,1000.];
plot()
colors = ["red","blue","Purple","brown"]
for (j,γ) in enumerate(γs)
    ud, J2, J3 = uucmols(γ)
    xx = simulate(ud)
    plot!(xx[1,:], xx[3,:], label="γ=$γ", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottom, color=colors[j])
end
scatter!(p[1,:],p[2,:], markershape=:star5, label="way-points", color="purple")
xlabel!("x₁")
ylabel!("x₃")
savefig("Final/p3f.png")
γs = 10.0.^(-1:0.1:150);
J2s = [];
J3s = [];
for (j,γ) in enumerate(γs)
    ud, J2, J3 = uucmols(γ)
    push!(J2s,J2)
    push!(J3s,J3)
end
plot(J2s, J3s, label="γ∈[10⁻¹,10¹⁵⁰]", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color = "green")
xlabel!("J₂")
ylabel!("J3")
savefig("Final/p3g.png")