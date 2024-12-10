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
include("simulator.jl")
print("\033c")

struct HistoryMonteCarloTreeSearch
    P # problem
    N # visit counts
    Q # action value estimates
    d # depth
    m # number of simulations
    c # exploration constant
    U # value function estimate
end

bonus(Nha, Nh) = Nha == 0 ? Inf : sqrt(log(Nh)/Nha)   # Algorithm 9.7 -> eqn. 22.3

function explore(π::HistoryMonteCarloTreeSearch, h)    # Algorithm 22.1
    A, N, Q, c = π.P.A, π.N, π.Q, π.c
    Nh = sum(get(N, (h,a), 0) for a in A)
    return argmax(a->Q[(h,a)] + c*bonus(N[(h,a)], Nh), A)
end

function simulate(π::HistoryMonteCarloTreeSearch, s, h, d)
    if d ≤ 0
        return π.U(s)
    end
    P, N, Q, c = π.P, π.N, π.Q, π.c
    S, A, TRO, γ = P.S, P.A, P.TRO, P.γ
    if !haskey(N, (h, first(A)))
        for a in A
            N[(h,a)] = 0
            Q[(h,a)] = 0.0
        end
        return π.U(s)
    end
    a = explore(π, h)
    s′, r, o = TRO(s,a)
    q = r + γ*simulate(π, s′, vcat(h, (a,o)), d-1)
    N[(h,a)] += 1
    Q[(h,a)] += (q-Q[(h,a)])/N[(h,a)]
    return q
end

function (π::HistoryMonteCarloTreeSearch)(b, h=[])
    for i in 1:π.m
        s = rand(SetCategorical(π.P.S, b))
        simulate(π, s, h, π.d)
    end
    return argmax(a->π.Q[(h,a)], π.P.A)
end


γ = 0.9;
S = 0:67227;
A = 1:6;
𝒪 = 0:67227;
P = POMDP(γ,S,A,𝒪,Transistion,Reward,ObsModel,TRO);
# N = Dict();
# Q = Dict();
d = 15;  # depth
m = 900;
c = 60;
# k_max = 10; # maximum number of iterations of QMDP  # Example 22.1
# πQMDP = solve(QMDP(k_max), P);
# save_object("QMPD.jld2",πQMDP)
# πQMDP = load_object("QMPD.jld2");
# U(b) = utility(πQMDP, b);
# αBAWS = baws_lowerbound(P);
# U(s) = αBAWS[s+1];
U(s) = 0.0
# HMCTS = HistoryMonteCarloTreeSearch(P,N,Q,d,m,c,U);

# b = zeros(length(S));
# T = 90;
# ss = Int.(zeros(T+1));
# ss[1] = state2id(0,[0,0,0,0,0]);
# b[ss[1]+1] = 1.0;
# h = [];
# for t in 1:T
#     a = HMCTS(b,h)
#     println(id2state(ss[t]))
#     println(a)
#     ss[t+1], r, o = TRO(ss[t],a)
#     h = vcat(h, (a,o))
#     if length(h)>3
#         deleteat!(h,1)
#     end
#     b = updateb(b,P,a,o)
#     if any(isnan.(b))
#         println(o)
#     end
# end

# xx = Int.(zeros(T+1));
# yy = Int.(zeros(T+1,5));
# for t in 1:T+1
#     xx[t],yy[t,:] = id2state(ss[t])
# end

T = 90;
K = 1;   # Number of trials
xx = Int.(zeros(T+1));
yy = Int.(zeros(T+1,5));
rr = zeros(T);
b = zeros(length(S));
b[ss[1]+1] = 1.0;
TT = zeros(K);
for k in 1:K
    t_init = time()
    N = Dict();
    Q = Dict();
    HMCTS = HistoryMonteCarloTreeSearch(P,N,Q,d,m,c,U);
    h = []
    ss = Int.(zeros(T+1))
    ss[:,1] .= state2id(0,[0,0,0,0,0])
    xtem = Int.(zeros(T+1))
    ytem = Int.(zeros(T+1,5))
    rtem = zeros(T)
    for t in 1:T
        a = HMCTS(b,h)
        println(id2state(ss[t]))
        println(a)
        ss[t+1], rtem[t], o = TRO(ss[t],a)
        xtem[t+1], ytem[t+1,:] = id2state(ss[t+1])
        h = vcat(h, (a,o))
        if length(h)>5
            deleteat!(h,1)
        end
        b = updateb(b,P,a,o)
        # b = zeros(length(S));
        # b[ss[t+1]+1] = 1.0;
    end
    rr += (rtem-rr)./k
    xx += (xtem-xx)./k
    yy += (ytem-yy)./k
    TT[k] = time()-t_init
    println("$(k)th trial done.")
end
println("Done, average time used = $(mean(TT)) sec")
println("Final accumulated reward = $(sum(rr))")
plotResults("HMCTS1",T,rr,xx,yy)