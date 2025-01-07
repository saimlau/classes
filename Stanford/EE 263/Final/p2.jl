# Problem 2
using LinearAlgebra
using Plots
print("\033c")

r = 10;
m = 15e3;
k = 15e4;
T = 30;
a = 0.1;
ω = 1.0;
f(t) = a*cos(ω*t);

K = diagm(1=>-ones(r-1)*k) + 2*k*I + diagm(-1=>-ones(r-1)*k);
K[r,r] = k;
G = zeros(r);
G[r] = 1;
A = [zeros(r,r) I;-K/m zeros(r,r)];
B = vcat(zeros(r),G/m);

h = 0.1;
Ad = exp(h*A);
Bd = inv(A)*(Ad-I)*B;
K = Int(T/h);
xd = zeros(2*r,K+1);
for k in 1:K
    xd[:,k+1] = Ad*xd[:,k]+Bd*f((k-1)*h)
end
plot((0:K).*h, xd[r,:], label="qᵣ", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottomright, color = "brown")
xlabel!("t (s)")
ylabel!("Horizontal Displacement (m)")
savefig("Final/p2g.png")