using CSV
using DataFrames
using DelimitedFiles
using LinearAlgebra
using Plots
using SparseArrays
using Random
using StatsBase
using Dierckx
using Distributions
using IterTools
using Base.Threads
using JLD2
include("utils.jl")
print("\033c")

γ = 0.9;
S = 0:67227;
A = 1:6;

T = 90;
K = 20   # Number of trials
xx = Int.(zeros(T+1));
yy = Int.(zeros(T+1,5));
# ss[1] = rand(S);
rr = zeros(T);
TT = zeros(K);
for k in 1:K
    t_init = time()
    ss = Int.(zeros(T+1));
    ss[:,1] .= state2id(0,[0,0,0,0,0]);
    xtem = Int.(zeros(T+1));
    ytem = Int.(zeros(T+1,5));
    rtem = zeros(T);
    for t in 1:T
        # x,y = id2state(ss[t])
        # a = x==0 ? 6 : argmin(y)
        a = rand(A)
        # println(id2state(ss[t]))
        # println(a)
        ss[t+1], rtem[t] = TR(ss[t],a)
        xtem[t+1], ytem[t+1,:] = id2state(ss[t+1])
    end
    rr += (rtem-rr)./k
    xx += (xtem-xx)./k
    yy += (ytem-yy)./k
    TT[k] = time()-t_init
    println("$(k)th trial done.")
end
println("Done, average time used = $(mean(TT)) sec")
println("Final accumulated reward = $(sum(rr))")
# plotResults("Baseline",T,rr,xx,yy)