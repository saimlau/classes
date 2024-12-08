using CSV
using DataFrames
using DelimitedFiles
using StatsPlots
using LinearAlgebra
using Plots
using SparseArrays
print("\033c")

function savePolicy(outfile, policy)
    writedlm(outfile,policy)
end

mutable struct QLearning
    S # state space (assumes 1:nstates)
    A # action space (assumes 1:nactions)
    γ # discount
    Q # action value function
    α # learning rate
end

function update!(model::QLearning, s, a, r, s′)
    γ, Q, α = model.γ, model.Q, model.α
    Q[s,a] += α*(r + γ*maximum(Q[s′,:]) - Q[s,a])
    return model
end

function interpolQ!(model::QLearning, ss, aa)
    for a in model.A
        s_idx = sort(unique(push!(ss[findall(aa.==a)],1,model.S[end])))
        for i in eachindex(s_idx)[1:end-1]
            n = s_idx[i+1]-s_idx[i]
            tem = Vector(0:n)./n
            model.Q[s_idx[i]:s_idx[i+1],a] = tem.*model.Q[s_idx[i+1],a].-(tem.-1).*model.Q[s_idx[i],a]
        end
    end
    return model
end

function compute(infile, iters)
    name = split(split(infile,"/")[end], ".")[1]
    outfile = name*".policy"
    data = CSV.read(infile, DataFrame)
    s = data[!,"s"]
    a = data[!,"a"]
    r = data[!,"r"]
    sp = data[!,"sp"]
    S = 1:302020
    A = 1:9
    Q = zeros(length(S), length(A)) # init value function
    α = 0.2 # learning rate
    model = QLearning(S,A,0.95,Q,α)
    last_policy = zeros(length(S))
    policy = ones(length(S))
    for i in 1:iters
        if i%100==0
            println("$i iterations done.")
        end
        # for t in sample(eachindex(s), length(s); replace=false)
        for t in eachindex(s)
            update!(model, s[t],a[t],r[t],sp[t])
        end
        interpolQ!(model, s, a)
        policy = [a[2] for a in argmax(model.Q,dims=2)]
        if norm(policy-last_policy) <= 1e-5
            println("Solved, i=$i.")
            break
        end
        last_policy = policy
    end
    savePolicy(outfile, policy)
end

t_init = time()
compute("data/large.csv", 6000)
t_f = time()
T = t_f-t_init
println("Done, time used = $T sec")
