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
N = Dict();
Q = Dict();
d = 15;  # depth
m = 50;
c = 1.2;
# k_max = 10; # maximum number of iterations of QMDP  # Example 22.1
# πQMDP = solve(QMDP(k_max), P);
# save_object("QMPD.jld2",πQMDP)
# πQMDP = load_object("QMPD.jld2");
# U(b) = utility(πQMDP, b);
αBAWS = baws_lowerbound(P);
U(s) = αBAWS[s+1];
HMCTS = HistoryMonteCarloTreeSearch(P,N,Q,d,m,c,U);

b = zeros(length(S));
T = 30;
ss = Int.(zeros(T+1));
ss[1] = state2id(3,[5,5,5,5,5]);
b[ss[1]] = 1.0;
for t in 1:T
    a = HMCTS(b)
    println(id2state(ss[t]))
    println(a)
    ss[t+1], r, o = TRO(ss[t],a)
    b = updateb(b,P,a,o)
    if any(isnan.(b))
        println(o)
    end
    # b = zeros(length(S));
    # b[ss[t+1]] = 1.0;
end