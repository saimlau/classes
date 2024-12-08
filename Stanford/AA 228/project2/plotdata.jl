using CSV
using DataFrames
using DelimitedFiles
using Plots
using StatsPlots
using LinearAlgebra
print("\033c")

data = CSV.read("data/small.csv", DataFrame)
s = data[!,"s"];
a = data[!,"a"];
r = data[!,"r"];
sp = data[!,"sp"];

xs = 1:10;
ys = 1:10;
S = 1:100;
A = 1:4;
directions = [([-0.5],[0.]),([0.5],[0.]),([0.],[0.5]),([0.],[-0.5])];
p = plot()
aa = [argmax([sum(a[s.==i].==j) for j in A]) for i in S]
for y in ys
    for x in xs
        scatter!([x],[y], label="", color="green")
        quiver!([x],[y],quiver=directions[aa[LinearIndices((10,10))[x,y]]], label="", color="red")
    end
end
xlims!(p,-0.5,10.5)
ylims!(p,-0.5,10.5)
savefig("small.png")

# medium
function sTopPosVel(s)
    Pos = (s-1)%500;
    Vel = Int((s-1-Pos)/500);
    return (Pos,Vel)
end

data = CSV.read("data/medium.csv", DataFrame);
s = data[!,"s"];
a = data[!,"a"];
r = data[!,"r"];
sp = data[!,"sp"];

S = 1:50000;
A = 1:7;
γ = 1.0;
Poses = 0:499;
Vels = 0:99;

ss = sTopPosVel.(s);
ss = transpose(reshape(reinterpret(Int, ss),2,:));
scatter(ss[:,1],ss[:,2])
savefig("medium_vel_vs_pos.png")
scatter(ss[:,1],r)
savefig("medium_r_vs_pos.png")

flag_pose = 0;
z_vel = 0;
j = 0;
idx = Vector{Int}();
for i in eachindex(s)
    if r[i]>5e4
        push!(idx,i)
        j += 1
        println(ss[i,:])
        flag_pose += (ss[i,1]-flag_pose)/j
        z_vel += (ss[i,2]-z_vel)/j
    end
    # if sp[i]!=s[i+1] && !in(i,[497, 747])
    #     println("$i")
    #     println("$(ss[i,:])")
    #     break
    # end
end
scatter(ss[idx,2],a[idx])
savefig("medium_a_vs_vel.png")
scatter(ss[2:end,2].-ss[1:end-1,2],a[1:end-1])
savefig("medium_a_vs_dv.png")
scatter(ss[2:end,1].-ss[1:end-1,1],ss[1:end-1,2])
savefig("medium_v_vs_dpos.png")

function calcR(s,a)
    ss = sTopPosVel.(s);
    ss = transpose(reshape(reinterpret(Int, ss),2,:))
    r = (100.0./abs.(ss[:,1].-460).^2)+(a-4).*(ss[:,1].-460)./3
    return r
end

R = stack([calcR(S,a) for a in A])
