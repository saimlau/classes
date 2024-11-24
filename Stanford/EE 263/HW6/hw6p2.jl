# Problem 2
using LinearAlgebra
using Plots
using LaTeXStrings

# Part a)
v(j) = [1 1 0 0;0 0 1 1]*[1 1 0 0;0 1 0 0;0 0 1 1;0 0 0 1]^j*[1/2 0;1 0;0 1/2;0 1];
A = zeros(6,142);
ts = [6,40,50];
for tt in 1:3
    for i in 2:ts[tt]
        A[tt*2-1:tt*2,(ts[tt]-i+1)*2-1:(ts[tt]-i+1)*2] = v(i-2)
    end
    A[tt*2-1:tt*2,ts[tt]*2-1:ts[tt]*2] = [0.5 0;0 0.5];
end
y_des = [1,-0.5,0,1,-1.5,0];
u_des = transpose(A)*inv(A*transpose(A))*y_des;
x = [0.]; y = [0.]; state = zeros(4);
for t in 1:70
    state = [1 1 0 0;0 1 0 0;0 0 1 1;0 0 0 1]*state+[0.5 0;1 0;0 0.5;0 1]*u_des[t*2-1:t*2]
    push!(x, state[1])
    push!(y, state[3])
end
fig = plot(x, y, label="Trajectory", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottom, color = "green", subplot=1, layout=(1,2), size = (660, 360))
scatter!(fig, [1,0,-1.5],[-0.5,1,0], label="way-points", color="red", markershape=:star5, subplot=1)
xlabel!(fig, "x", subplot=1)
ylabel!(fig, "y", subplot=1)
plot!(fig, 0:70, u_des[1:2:end], subplot=2, label="u₁", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright)
plot!(fig, 0:70, u_des[2:2:end], subplot=2, label="u₂", linewidth=2)
xlabel!(fig, "t", subplot=2)
savefig("HW6/p2a.png")

# Part b)
J1 = [];
J2 = [];
mu = 10.0.^(-1:0.01:5);
for m in mu
    u = inv(transpose(A)*A+m*I)*transpose(A)*y_des
    push!(J1, norm(A*u-y_des)^2)
    push!(J2, norm(u)^2)
end
plot(J2,J1, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color = "purple", label="μ∈[10^-1,10^5]")
xlabel!("J1")
ylabel!("J2")
savefig("HW6/p2b_tradeoff.png");
scatter([1,0,-1.5],[-0.5,1,0], label="way-points", color="red", markershape=:star5, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottom)
xlabel!("x")
ylabel!("y")
for p in -2:2:10
    m = 10^(p/2)
    u = inv(transpose(A)*A+m*I)*transpose(A)*y_des
    x = [0.]; y = [0.]; state = zeros(4);
    for t in 1:70
        state = [1 1 0 0;0 1 0 0;0 0 1 1;0 0 0 1]*state+[0.5 0;1 0;0 0.5;0 1]*u[t*2-1:t*2]
        push!(x, state[1])
        push!(y, state[3])
    end
    plot!(x, y, label="μ=10^($p/2)", linewidth=2)
end
savefig("HW6/p2c.png")

# Part g
function init_W()
    W = zeros(142,142);
    for t in 0:70
        tt = t+1
        for i in 2:t
            if t>=2
                W[tt*2-1:tt*2,(t-i+1)*2-1:(t-i+1)*2] = v(i-2)
            end
        end
        if t>=1
            W[tt*2-1:tt*2,t*2-1:t*2] = [0.5 0;0 0.5]
        end
    end
    return W
end
scatter([1,0,-1.5],[-0.5,1,0], label="way-points", color="red", markershape=:star5, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottomleft)
xlabel!("x")
ylabel!("y")
for p in 0:2:20
    γ = 10^(p/2)
    W = init_W()+sqrt(γ)*I
    Σ = transpose(W)*W;
    u_opt = inv(Σ)*transpose(A)*inv(A*inv(Σ)*transpose(A))*y_des
    qt = init_W()*u_opt
    plot!(qt[1:2:end], qt[2:2:end], label="γ=10^($p/2)", linewidth=2)
end
scatter!([1,0,-1.5],[-0.5,1,0], label="", color="red", markershape=:star5, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottomleft)
xlims!(-6,3)
ylims!(-3,1.5)
savefig("HW6/p2g.png")

# Part h
function circlesShape(r)
    θ = LinRange(0, 2*pi, 500)
    r*sin.(θ), r*cos.(θ)
end
plot(circlesShape(2), seriestype=[:shape],lw=0.5,c=RGBA{Float64}(0, 100, 100, 0.3),linecolor=:black, label="radio range")
scatter!([1,0,-1.5],[-0.5,1,0], label="way-points", color="red", markershape=:star5, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottomleft)
xlabel!("x")
ylabel!("y")
p = 6
γ = 10^(p/2)
W = init_W()+sqrt(γ)*I
Σ = transpose(W)*W;
u_opt = inv(Σ)*transpose(A)*inv(A*inv(Σ)*transpose(A))*y_des
qt = init_W()*u_opt
plot!(qt[1:2:end], qt[2:2:end], label="γ=10^($p/2)", linewidth=2)
savefig("HW6/p2h.png")

# Part i
J2 = [];
J3 = [];
γs = 10.0.^((6:0.1:20)./2);
for γ in γs
    W = init_W()+sqrt(γ)*I
    Σ = transpose(W)*W;
    u = inv(Σ)*transpose(A)*inv(A*inv(Σ)*transpose(A))*y_des
    qt = init_W()*u
    push!(J2, norm(u)^2)
    push!(J3, norm(qt)^2)
end
plot(J3,J2, framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
linewidth=2, background_color_legend = nothing, foreground_color_legend = nothing, legend=:topright, color = "purple", label="γ∈[10^3,10^10]")
xlabel!("J2")
ylabel!("J3")
savefig("HW6/p2i_tradeoff.png")