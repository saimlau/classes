# Problem 5
using LinearAlgebra
using Plots
using LaTeXStrings
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW6/nba_ranking_data.json");
records = data["records"]
scores = data["scores"]
T = data["T"]
teams = data["teams"]
n = data["n"]

A = zeros(size(scores));
idx = scores.>0;
A[idx] = scores[idx];
z0 = records[:,1]./(sum(records,dims=2));
x0 = z0/norm(z0);
L, V = eigen(A);
v1 = V[:,argmax(abs.(L))];
L1 = L[argmax(abs.(L))];
L2 = L[L.!=L1][argmax(abs.(L[L.!=L1]))];
L3 = L[L.!=L1 .&& L.!=L2][argmax(abs.(L[L.!=L1 .&& L.!=L2]))];
x_bar = sign(dot(v1,x0))*v1;
e = [norm(x0-x_bar)];
X = reshape(x0,:,1);
for t=1:T
    X = hcat(X,A*X[:,end]./norm(A*X[:,end]))
    push!(e, norm(X[:,end]-x_bar))
end
plot(0:T,X[10,:], framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, label="j=10")
plot!(0:T, X[20,:], lw=2, label="j=20")
plot!(0:T, X[30,:], lw=2, label="j=30")
scatter!([T],[real(x_bar[10])], label="x̄₁₀", color="red", markershape=:star5, ms=9)
scatter!([T],[real(x_bar[20])], label="x̄₂₀", color="yellow", markershape=:star5, ms=9)
scatter!([T],[real(x_bar[30])], label="x̄₃₀", color="green", markershape=:star5, ms=9)
xlabel!("t")
savefig("HW6/p5ci.png")

ρ = abs(L2/L1);
plot(0:T, e, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, label="||e(t)||", yaxis=:log)
xlabel!("t")
plot!(1:T,ρ.^(1:T), label="slope=log₁₀(ρ)=$(round(log(10,ρ), sigdigits=5))", ls = :dash, lw=2)
savefig("HW6/p5cii.png")

overrated = teams[argmax(real.(x0-x_bar))];
underrated = teams[argmax(real.(x_bar-x0))];
println("Would be most overrated: $overrated")
println("Would be most underrated: $underrated")
