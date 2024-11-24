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
include("utils.jl")
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

function explore(π::HistoryMonteCarloTreeSearch, h)
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

function compute()
    γ = 0.9
    S = 
    P = POMDP()
end