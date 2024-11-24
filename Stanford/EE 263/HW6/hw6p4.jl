# Problem 4
using LinearAlgebra
using Plots
using LaTeXStrings
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW6/power_iteration_convergence_rate_data.json");
A = data["A"]
T = data["T"]
x0 = data["x0"]

L, V = eigen(A);
v1 = V[:,3];
x_bar = sign(dot(v1,x0))*v1;
e = [norm(x0-x_bar)];
X = reshape(x0,(:,1));
for t in 1:T
    X = hcat(X,A*X[:,end]./norm(A*X[:,end]))
    push!(e, norm(X[:,end]-x_bar))
end
plot(0:T,X[1,:], framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottom, label="x₁")
plot!(0:T, X[2,:], lw=2, label="x₂")
plot!(0:T, X[3,:], lw=2, label="x₃")
scatter!([T],[x_bar[1]], label="x̄₁", color="red", markershape=:star5, ms=9)
scatter!([T],[x_bar[2]], label="x̄₂", color="yellow", markershape=:star5, ms=9)
scatter!([T],[x_bar[3]], label="x̄₃", color="green", markershape=:star5, ms=9)
xlabel!("t")
savefig("HW6/p4c.png")

plot(0:T, e, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, label="||e(t)||", yaxis=:log)
xlabel!("t")
plot!(1:T,abs(L[1]/L[3]).^(1:T), label=L"slope=\log_{10}\left\vert\frac{\lambda_2}{\lambda_1}\right\vert", ls = :dash, lw=2)
savefig("HW6/p4d.png");