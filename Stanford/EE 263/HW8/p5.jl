# Problem 5
using LinearAlgebra
using Plots
print("\033c")
A = [0.5 0.7 -0.9 -0.5;0.4 -0.7 0.1 0.3;0.7 0 -0.6 0.1;0.4 -0.1 0.8 -0.5];
B = [1,1,0,0];
n = length(B);
x_des = [0.8,2.3,-0.7,-0.3];

function CC(T)
    CT = zeros(n,0)
    for t in 0:T-1
        CT = hcat(CT,A^t*B)
    end
    return CT
end

for T in 1:30
    CT = CC(T)
    u_est = pinv(CT)*x_des
    if norm(CT*u_est-x_des)<1e-5
        println("Smallest T = $T")
        display(CT)
        display(u_est) 
        break
    end
end
rank(CC(4))
Tmin = 4;
E = [];
x_des = [-1,1,0,1];
for T in Tmin:30
    CT = CC(T)
    u = transpose(CT)*inv(CT*transpose(CT))*x_des
    push!(E,norm(u)^2)
    # println(round.(CT*u,sigdigits=3))
end
plot(Tmin:30, E, label="Eₘᵢₙ(T)", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color = "brown")
xlabel!("T")
ylabel!("E(T)")
savefig("HW8/p5c.png");