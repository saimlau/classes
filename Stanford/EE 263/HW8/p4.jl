# Problem 4
using LinearAlgebra
using Symbolics
using Plots
print("\033c")

# Scheme 1
tem = diagm(1=>ones(4));
tem[diagind(tem)] = -ones(5);
A = [zeros(5,5) I;tem -I];
b = [0,0,0,0,0,-1,-1,-1,-1,5];
L1, ~ = eigen(A);
tem2 = inv(A)*b;
println("Scheme 1:")
println("\t R{λ(A)}=$(round.(real.(L1),sigdigits=2))")
println("\t A⁻¹b=$tem2")

# Scheme 2
tem = diagm(1=>ones(4)/2);
tem[diagind(tem)] = -ones(5);
tem += diagm(-1=>[.5,.5,.5,0]);
tem[1,2] = 1;
B = [zeros(5,5) I;tem -I];
d = [0,0,0,0,0,-1,0,0,0,5];
L2, ~ = eigen(B);
tem2 = inv(B)*d;
println("Scheme 2:")
println("\t R{λ(B)}=$(round.(real.(L2),sigdigits=2))")
println("\t B⁻¹d=$tem2")

# Scheme 3
C = [zeros(5,5) I;-I -I];
h = [0,0,0,0,0,1,2,3,4,5];
L3, ~ = eigen(C);
tem2 = inv(C)*h;
println("Scheme 3:")
println("\t R{λ(C)}=$(round.(real.(L3),sigdigits=2))")
println("\t C⁻¹h=$tem2")

# Collisions
function find_collision(ss,ts)
    tem = (zeros(length(ss[:,1])))
    for i in eachindex(tem)
        tem2 = findfirst(ss[i,:].<=0)
        if !isnothing(tem2)
            tem[i] = abs(ss[i,tem2])<ss[i,tem2-1] ? tem2 : tem2-1
        else
            tem[i] = Inf
        end
    end
    tem3 = argmin(tem)
    try
        tem4 = ts[Int(minimum(tem))]
        println("$(tem3+1) and $(tem3) collide at t ≈ $tem4")
    catch
        tem4 = argmin(ss)
        println("$(tem4[1]+1) and $(tem4[1]) reached closest dis. $(round(ss[tem4], sigdigits=5)) at t ≈ $(ts[tem4[2]])")
    end
end
colors = ["red","blue","purple","green"]
s(A,t) = [-1 1 0 0 0;0 -1 1 0 0;0 0 -1 1 0;0 0 0 -1 1]*((exp(A*t)[1:5,1:5])*[-1,0,0,1,2]+[1,2,3,4,5]);
ts = 0:0.01:10;
ss = ones(4,length(ts));
for i in eachindex(ts)
    ss[:,i] = s(A,ts[i])
end
plot()
for i in 1:4
    plot!(ts,ss[i,:], linewidth=2, color=colors[i], label="s$i")
end
plot!()
title!("Scheme 1")
savefig("HW8/p4b1.png");
find_collision(ss, ts);

ss = ones(4,length(ts));
for i in eachindex(ts)
    ss[:,i] = s(B,ts[i])
end
plot()
for i in 1:4
    plot!(ts,ss[i,:], linewidth=2, color=colors[i], label="s$i")
end
plot!()
title!("Scheme 2")
savefig("HW8/p4b2.png");
find_collision(ss, ts);

ss = ones(4,length(ts));
for i in eachindex(ts)
    ss[:,i] = s(C,ts[i])
end
plot()
for i in 1:4
    plot!(ts,ss[i,:], linewidth=2, color=colors[i], label="s$i")
end
plot!()
title!("Scheme3")
savefig("HW8/p4b3.png");
find_collision(ss, ts);