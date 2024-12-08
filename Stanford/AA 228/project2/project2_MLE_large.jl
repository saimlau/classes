using CSV
using DataFrames
using DelimitedFiles
using StatsPlots
using LinearAlgebra
using Plots
using SparseArrays
using Random
using StatsBase
using Dierckx
print("\033c")

global updated_states = []

struct MDP   # Algorithm 7.1
    γ # discount factor
    S # state space
    A # action space
    T # transition function
    R # reward function
    # TR # sample transition and reward
    planner
end

mutable struct MaximumLikelihoodMDP   # Algorithm 16.1
    S # state space (assumes 1:nstates)
    A # action space (assumes 1:nactions)
    N # transition count N(s,a,s′)
    ρ # reward sum ρ(s, a)
    γ # discount
    U # value function
    planner
end

struct ValueFunctionPolicy   # Algorithm 7.5
    P # problem
    U # utility function
end

struct ValueIteration   # Algorithm 7.8
    k_max # maximum number of iterations
end;

struct RandomizedUpdate  # Algorithm 16.4
    m # number of updates
end

function lookahead(P::MDP, U, s, a)
    S, T, R, γ = P.S, P.T, P.R, P.γ
    tem = findall(N.!=0)
    tem2 = copy(N)
    # print(size(tem)) # Gave 4498
    for i in sample(tem, P.planner.m; replace=false)
        indx = unique(idxss[:,3][findall(eachrow(idxss[:,1:2]).==[[i[1],i[2]]])])
        tem2[i[1],i[2]] = dot(T.(i[1],i[2],indx), U[indx])
        push!(updated_states, i[1])
    end
    unique!(updated_states)
    return R + γ*tem2
end

function backup(P::MDP, U, s)    # Algoritm 7.7
    return maximum(lookahead(P, U, s, P.A), dims=2)
end

function greedy(P::MDP, U, s)  # Algorithm 7.5
    return reshape(argmax(lookahead(P, U, s, P.A), dims=2), length(P.S))
end
(π::ValueFunctionPolicy)(S) = greedy(π.P, π.U, S)

function interpolU(U, S)
    fd = Spline1D(updated_states, U[updated_states], bc="nearest")
    return fd(S)
end

function solve(M::ValueIteration, P::MDP)   # Algorithm 7.8
    U = [0.0 for s in P.S]
    for k = 1:M.k_max
        if k%10==0
            println("$k iterations on U done.")
        end
        global updated_states = []
        U = reshape(backup(P, U, P.S), length(P.S))
        push!(updated_states, 1, P.S[end])
        sort!(unique!(updated_states))
        U = interpolU(U, P.S)
    end
    return ValueFunctionPolicy(P, U)
end;

function MDP(model::MaximumLikelihoodMDP)   # Algorithm 16.2
    N, ρ, S, A, γ = model.N, model.ρ, model.S, model.A, model.γ
    T(s,a,s′) = N[s,a]==0 ? 0.0 : N_fn(s,a,s′)/N[s,a];
    R = copy(ρ)./N
    R[isnan.(R)] .= 0.0
    return MDP(γ,S,A,T,R, model.planner)
end;

function policy_evaluation(P::MDP, π)  # Algorithm 7.4
    S, R, T, γ = P.S, P.R, P.T, P.γ
    R′ = [R(s, π(s)) for s in S]
    T′ = [T(s, π(s), s′) for s in S, s′ in S]
    return (I - γ*T′)\R′
end;

function savePolicy(outfile, policy)
    writedlm(outfile,policy)
end

function update!(M::ValueIteration, planner::RandomizedUpdate, model)   # Algorithm 16.3
    P = MDP(model)
    println("Now solving for U.")
    sol = solve(M, P)
    U = sol.U
    copy!(model.U, U)
    return sol
end;

t_init = time()
data = CSV.read("data/large.csv", DataFrame);
s = data[!,"s"];
a = data[!,"a"];
r = data[!,"r"];
sp = data[!,"sp"];

S = 1:302020
A = 1:9
γ = 0.95;
planner = RandomizedUpdate(3000);
U = zeros(length(S));
tem_data = hcat(s,a,sp);
global idxss = unique(tem_data,dims=1);
N_counts = [sum(eachrow(tem_data).==[ℓ]) for ℓ in eachrow(idxss)];
N_fn(s,a,s′) = in([s,a,s′],eachrow(idxss)) ? N_counts[eachrow(idxss).==[[s,a,s′]]][1] : 0;
N = zeros(length(S),length(A));
idxsss = unique(idxss[:,1:2], dims=1);
for i in eachrow(idxsss)
    N[i[1],i[2]] = sum(eachrow(tem_data[:,1:2]).==[i])
end

idxs = unique(hcat(s,a),dims=1);
ρ = sparse(idxs[:,1],idxs[:,2],zeros(length(idxs[:,1])),length(S),length(A));
for i in eachindex(s)
    ρ[s[i],a[i]] += r[i];
end

model = MaximumLikelihoodMDP(S,A,N,ρ,γ,U,planner);
M = ValueIteration(600)
π = update!(M, model.planner, model)
policy = π(S);
policy = [p[2] for p in policy];
outfile = "large.policy"
savePolicy(outfile, policy)
t_f = time()
T = t_f-t_init
println("Done, time used = $T sec")

