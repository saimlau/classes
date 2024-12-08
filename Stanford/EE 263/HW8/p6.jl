# Problem 6
using LinearAlgebra
using Plots
print("\033c")
P = [-1 -.5 0 0.5 1 0; 2 0 1 0 2 0]; s = [3; 7; 10; 13; 17; 20];
m = length(P[1,:]);

A = diagm(2=>[1,1]);
B = [0 0;0 0;1 0;0 1];
C4 = [B A*B A^2*B A^3*B];

# Part d
h = 0.1;
T = s[m];
n = Int(T/h);
Ad = I+A*h;
Bd = [h^2/2 0;0 h^2/2;1 0;0 1];
Cd = I(4);
vP = vec(P);
function CC(K)
    CK = zeros(4,0)
    for k in 0:K-1
        CK = hcat(CK,Ad^k*Bd)
    end
    return CK
end
C = zeros(2*m,2*Int(n));
for i in 1:m
    C[2*i-1:2*i,:] = [zeros(2,Int(2*(n-s[i]/h))) CC(Int(s[i]/h))[1:2,:]]
end
ud = transpose(C)*inv(C*transpose(C))*vP;
ud = reshape(ud,2,:);
ud = reverse(ud,dims=2);
J1 = norm(ud)^2;
y = zeros(4,n+1);
for k in 1:Int(n)
    y[:,k+1] = Ad*y[:,k]+Bd*ud[:,k]
end
plot(framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:top)
plot!(y[1,:],y[2,:], label="trajectory", linewidth=2, color="brown")
scatter!(P[1,:],P[2,:], markerstyle="star5", label="P", color="purple")
xlabel!("q₁")
ylabel!("q2")
savefig("HW8/p6d_traj.png")
plot(framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright)
plot!((0:n-1)*h, ud[1,:], label="ud₁", linewidth=2, color="red")
plot!((0:n-1)*h, ud[2,:], label="ud₁", linewidth=2, color="blue")
xlabel!("t")
ylabel!("")
title!("J₁≈$(round(J1,sigdigits=3))")
savefig("HW8/p6d_ctrl.png")

# Part e
C2 = zeros(4*m,2*Int(n));
for i in 1:m
    C2[4*i-3:4*i,:] = [zeros(4,Int(2*(n-s[i]/h))) CC(Int(s[i]/h))]
end
vP2 = [P[:,1];0;0];
for i in 2:m
    vP2 = vcat(vP2,[P[:,i];0;0])
end
ud = transpose(C2)*inv(C2*transpose(C2))*vP2;
ud = reshape(ud,2,:);
ud = reverse(ud,dims=2);
J1 = norm(ud)^2;
y = zeros(4,n+1);
for k in 1:Int(n)
    y[:,k+1] = Ad*y[:,k]+Bd*ud[:,k]
end
plot(framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:top)
plot!(y[1,:],y[2,:], label="trajectory", linewidth=2, color="brown")
scatter!(P[1,:],P[2,:], markerstyle="star5", label="P", color="purple")
xlabel!("q₁")
ylabel!("q2")
savefig("HW8/p6e_traj.png")
plot(framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright)
plot!((0:n-1)*h, ud[1,:], label="ud₁", linewidth=2, color="red")
plot!((0:n-1)*h, ud[2,:], label="ud₁", linewidth=2, color="blue")
xlabel!("t")
ylabel!("")
title!("J₁≈$(round(J1,sigdigits=3))")
savefig("HW8/p6e_ctrl.png")

# Part f
μ = 100;
F = (I(2*n)+diagm(2=>-ones(2*n-2)))[1:end-2,:];
R = [I;sqrt(μ)*F];
RR = [transpose(R)*R transpose(C2);C2 zeros(24,24)];
udL = pinv(RR)*[zeros(2*n);vP2];
ud = udL[1:2*n];
ud = reshape(ud,2,:);
ud = reverse(ud,dims=2);
J1 = norm(ud)^2;
J2 = norm(F*udL[1:2*n])^2;
y = zeros(4,n+1);
for k in 1:Int(n)
    y[:,k+1] = Ad*y[:,k]+Bd*ud[:,k]
end
plot(framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:top)
plot!(y[1,:],y[2,:], label="trajectory", linewidth=2, color="brown")
scatter!(P[1,:],P[2,:], markerstyle="star5", label="P", color="purple")
xlabel!("q₁")
ylabel!("q2")
savefig("HW8/p6f_traj.png")
plot(framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:top)
plot!((0:n-1)*h, ud[1,:], label="ud₁", linewidth=2, color="red")
plot!((0:n-1)*h, ud[2,:], label="ud₁", linewidth=2, color="blue")
xlabel!("t")
ylabel!("")
title!("μ=$μ, J₁≈$(round(J1,sigdigits=3)), J₂≈$(round(J2,sigdigits=3))")
savefig("HW8/p6f_ctrl.png")

J1s = [];
J2s = [];
for μ in 0:1:300
    F = (I(2*n)+diagm(2=>-ones(2*n-2)))[1:end-2,:];
    R = [I;sqrt(μ)*F];
    RR = [transpose(R)*R transpose(C2);C2 zeros(24,24)];
    udL = pinv(RR)*[zeros(2*n);vP2];
    ud = udL[1:2*n];
    ud = reshape(ud,2,:);
    ud = reverse(ud,dims=2);
    push!(J1s, norm(ud)^2);
    push!(J2s, norm(F*udL[1:2*n])^2);
end
plot(J2s, J1s, label="μ∈[0,300]", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
        linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color="green")
xlabel!("J2")
ylabel!("J1")
savefig("HW8/p6f_trade-off.png")