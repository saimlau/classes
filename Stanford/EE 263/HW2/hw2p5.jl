# problem 5
using LinearAlgebra

n = 20;
ti = 0.96.*ones(n-1);
ri = 0.02.*ones(n-1);

function find_S(n, ti, ri)
    tem1 = diagm(-1 => ti);
    tem1[1,1] = 1;
    tem2 = diagm(0 => pushfirst!(copy(ri),0));
    tem3 = diagm(0 => push!(copy(ri),1));
    tem4 = diagm(1 => ti);

    A = vcat([tem1 tem2; tem3 tem4]);

    e_val,e_vec = eigen(A);
    indx = findfirst(x -> x==1,e_val)
    z = e_vec[:,indx];
    S = z[n+1]/z[1];
    return S
end

S = find_S(n, ti, ri)
println("S = $S")
Ss = [];

for k=1:(n-1)
    ti_tem = copy(ti);
    ri_tem = copy(ri);
    ti_tem[k] = 0.02;
    ri_tem[k] = 0.96;
    push!(Ss, find_S(n,ti_tem,ri_tem))
end

k_sol = 0;
min_k = 999.0;
for k=1:(n-1)
    if abs(Ss[k]-0.7)<min_k 
        k_sol=k;
        min_k=abs(Ss[k]-0.7)
    end
end
println("k = $k_sol")

