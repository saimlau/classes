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

struct MDP           # Algorithm 7.1
    γ # discount factor
    S # state space
    A # action space
    T # transition function
    R # reward function
    TR # sample transition and reward
end

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

function TR(s,a)                        # x in {0,1,2,3}
    x,y = id2state(s)
    r = Reward(s,a)
    if a ==6 #|| (x==0 && rand(SetCategorical([0,1],[0.7,0.3]))==1)
        x′ = rand([1,2,3])       
        y′ = [yi==6 ? rand(dis2) : yi-1 for yi in y]
    else
        r -= 3
        if x==0
            x′ = x
            y′ = copy(y)
            y′[a] = rand(SetCategorical([6,y[a]],[0.2,0.8]))
        else
            x′ = x+rand(dis3)
            y′ = copy(y)
            y′[a] = 6
        end
    end
    y′[y′.<0] .= 0
    # if !(state2id(x′,y′) in possTransit(s,a))
    #     println(x′)
    #     println(y′)
    #     error("Getting impossible state!!!!!")
    # end
    return (state2id(x′,y′), r)
end

γ = 0.9;
S = 0:67227;
A = 1:6;
P = MDP(γ,S,A,Transistion,Reward,TR)
N = Dict();
Q = Dict();
d = 30;  # depth
m = 100;
c = 1.2;
U(s) = 0.0
MCTS = MonteCarloTreeSearch(P,N,Q,d,m,c,U)


T = 300;
ss = Int.(zeros(T+1));
ss[1] = state2id(3,[0,0,0,0,0]);
rr = zeros(T)
for t in 1:T
    a = MCTS(ss[t])
    println(id2state(ss[t]))
    println(a)
    ss[t+1], rr[t] = TR(ss[t],a)
end
plot(1:T,rr)
