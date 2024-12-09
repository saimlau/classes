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

struct MonteCarloTreeSearch
    P # problem
    N # visit counts
    Q # action value estimates
    d # depth
    m # number of simulations
    c # exploration constant
    U # value function estimate
end
function (π::MonteCarloTreeSearch)(s)
    for k in 1:π.m
        simulate!(π, s)
    end
    return argmax(a->π.Q[(s,a)], π.P.A)
end

function simulate!(π::MonteCarloTreeSearch, s, d=π.d)
    if d ≤ 0
        return π.U(s)
    end
    P, N, Q, c = π.P, π.N, π.Q, π.c
    A, TR, γ = P.A, P.TR, P.γ
    if !haskey(N, (s, first(A)))
        for a in A
            N[(s,a)] = 0
            Q[(s,a)] = 0.0
        end
        return π.U(s)
    end
    a = explore(π, s)
    s′, r = TR(s,a)
    q = r + γ*simulate!(π, s′, d-1)
    N[(s,a)] += 1
    Q[(s,a)] += (q-Q[(s,a)])/N[(s,a)]
    return q
end

bonus(Nsa, Ns) = Nsa == 0 ? Inf : sqrt(log(Ns)/Nsa)

function explore(π::MonteCarloTreeSearch, s)
    A, N, Q, c = π.P.A, π.N, π.Q, π.c
    Ns = sum(N[(s,a)] for a in A)
    return argmax(a->Q[(s,a)] + c*bonus(N[(s,a)], Ns), A)
end

γ = 0.9;
S = 0:67227;
A = 1:6;
P = MDP(γ,S,A,Transistion,Reward,TR)
d = 30;  # depth
m = 150;
c = 1.2;
U(s) = 0.0


T = 90;
K = 20   # Number of trials
xx = Int.(zeros(T+1));
yy = Int.(zeros(T+1,5));
# ss[1] = rand(S);
rr = zeros(T);
TT = zeros(K);
for k in 1:K
    t_init = time()
    N = Dict();
    Q = Dict();
    MCTS = MonteCarloTreeSearch(P,N,Q,d,m,c,U)
    ss = Int.(zeros(T+1));
    ss[:,1] .= state2id(0,[0,0,0,0,0]);
    xtem = Int.(zeros(T+1));
    ytem = Int.(zeros(T+1,5));
    rtem = zeros(T);
    for t in 1:T
        a = MCTS(ss[t])
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
plotResults("MCTS",T,rr,xx,yy)

