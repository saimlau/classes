using CSV
using DataFrames
using DelimitedFiles
using StatsPlots
using LinearAlgebra
using Plots
using SparseArrays
print("\033c")

struct MDP
    γ # discount factor
    S # state space
    A # action space
    T # transition function
    R # reward function
    # TR # sample transition and reward
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

struct FullUpdate end   # Algorithm 16.3

struct ValueIteration   # Algorithm 7.8
    k_max # maximum number of iterations
end;

function lookahead(P::MDP, U, s, a)
    S, T, R, γ = P.S, P.T, P.R, P.γ
    return R[s,a] + γ*sum(T[s][a,s′]*U[s′] for s′ in S)
end

function backup(P::MDP, U, s)    # Algoritm 7.7
    return maximum(lookahead(P, U, s, a) for a in P.A)
end

function greedy(P::MDP, U, s)  # Algorithm 7.5
    u, a = findmax(a->lookahead(P, U, s, a), P.A)
    return (a=a, u=u)
end
(π::ValueFunctionPolicy)(s) = greedy(π.P, π.U, s).a

function solve(M::ValueIteration, P::MDP)   # Algorithm 7.8
    U = [0.0 for s in P.S]
    for k = 1:M.k_max
        U = [backup(P, U, s) for s in P.S]
    end
    return ValueFunctionPolicy(P, U)
end;

function MDP(model::MaximumLikelihoodMDP)   # Algorithm 16.2
    N, ρ, S, A, γ = model.N, model.ρ, model.S, model.A, model.γ
    T, R = copy(N), copy(ρ)
    for s in S
        for a in A
            n = sum(N[s][a,:])
            if n == 0
                T[s][a,:] .= 0.0
                R[s,a] = 0.0
            else
                T[s][a,:] = N[s][a,:] / n
                R[s,a] = ρ[s,a] / n
            end
        end
    end
    return MDP(γ,S,A,T,R)
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

function update!(M::ValueIteration, planner::FullUpdate, model)   # Algorithm 16.3
    P = MDP(model)
    sol = solve(M, P)
    U = sol.U
    copy!(model.U, U)
    return sol
end;


t_init = time()
data = CSV.read("data/small.csv", DataFrame);
s = data[!,"s"];
a = data[!,"a"];
r = data[!,"r"];
sp = data[!,"sp"];

S = 1:100;
A = 1:4;
γ = 0.95;
planner = FullUpdate();
U = zeros(length(S))
N = [];
for i in S
    idxs = unique(hcat(a[s.==i],sp[s.==i]),dims=1)
    Ni = sparse(idxs[:,1],idxs[:,2],zeros(length(idxs[:,1])),length(A),length(S))
    for j in findall(s.==i)
        Ni += sparse([a[j]],[sp[j]],[1],length(A),length(S))
    end
    push!(N,Ni)
end
idxs = unique(hcat(s,a),dims=1);
ρ = sparse(idxs[:,1],idxs[:,2],zeros(length(idxs[:,1])),length(S),length(A));
for i in eachindex(s)
    ρ[s[i],a[i]] += r[i];
end

model = MaximumLikelihoodMDP(S,A,N,ρ,γ,U,planner);
M = ValueIteration(300)
π = update!(M, model.planner, model)
policy = [π(i) for i in S]
outfile = "small.policy"
savePolicy(outfile, policy)
t_f = time()
T = t_f-t_init
println("Done, time used = $T sec")
