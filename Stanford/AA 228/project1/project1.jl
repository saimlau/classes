using Graphs
using Printf
using DataFrames
using CSV
using IterTools
using Random
using SpecialFunctions
using LinearAlgebra
using GraphPlot
using Compose, Cairo, Fontconfig

"""
    write_gph(dag::DiGraph, idx2names, filename)

Takes a DiGraph, a Dict of index to names and a output filename to write the graph in `gph` format.
"""
function write_gph(dag::DiGraph, idx2names, filename)
    open(filename, "w") do io
        for edge in edges(dag)
            @printf(io, "%s,%s\n", idx2names[src(edge)], idx2names[dst(edge)])
        end
    end
    node_names = [idx2names[i] for i in 1:nv(dag)]
    p = gplot(dag; nodelabel=node_names)
    fname = rsplit(filename, ".")[1]
    draw(PDF(fname*".pdf", 16cm, 16cm), p)
end

function prior(n::Int, r::Vector{Int64}, G::DiGraph)
    q = [prod([r[j] for j in inneighbors(G,i)]) for i in 1:n]
    return [ones(q[i], r[i]) for i in 1:n]
end

function init_M(A)
    return [zeros(size(a)) for a in A]
end

function statistics_update(old_M, )
    return
end

function sub2ind(siz, x)
    k = vcat(1, cumprod(siz[1:end-1]))
    return dot(k, x .- 1) + 1
end

function statistics(n, r, G, D::DataFrame, M_old, i_affected)
    q_i = prod([r[j] for j in inneighbors(G,i_affected)])
    # M = init_M(prior(n, r, G))
    # M = [zeros(size(a)) for a in M_old]
    # M = [M_old[i] for i in 1:n]
    # M[i_affected] = zeros(q_i, r[i_affected])
    # for o in eachrow(D)
    #     k = o[i_affected]
    #     parents = inneighbors(G,i_affected)
    #     j = 1
    #     if !isempty(parents)
    #         j = sub2ind(r[parents], Vector(o[parents]))
    #     end
    #     M[i_affected][j,k] += 1.0
    # end

    M = init_M(prior(n, r, G))
    for o in eachrow(D)
        for i in 1:n
            k = o[i]
            parents = inneighbors(G,i)
            j = 1
            if !isempty(parents)
                j = sub2ind(r[parents], Vector(o[parents]))
            end
            M[i][j,k] += 1.0
        end
    end
    return M
end

function bayesian_score_component(M, A)
    p = sum(loggamma.(A + M))
    p -= sum(loggamma.(A))
    p += sum(loggamma.(sum(A,dims=2)))
    p -= sum(loggamma.(sum(A,dims=2) + sum(M,dims=2)))
    return p
end

function bayesian_score(n::Int, r::Vector{Int64}, G::DiGraph, D::DataFrame, M, i_affected)
    M = statistics(n, r, G, D, M, i_affected)
    alpha = prior(n, r, G)
    return sum(bayesian_score_component(M[i], alpha[i]) for i in 1:n)
end

function compute(infile, outfile)
    dag = DiGraph()
    idx2names = Dict()
    data = CSV.read(infile, DataFrame)
    n = length(names(data))
    add_vertices!(dag, n)
    for i=1:n
        idx2names[i] = names(data)[i]
    end
    # r = [length(unique(data[:,i])) for i in 1:n]
    r = [maximum(data[:,i]) for i in 1:n]
    # M = init_M(prior(n, r, dag))
    
    # Max_itr = 1000
    ordering = shuffle(shuffle(1:n))
    M = init_M(prior(n, r, dag))

    for (k,i) in enumerate(ordering[2:end])
        y = bayesian_score(n, r, dag, data, M, ordering[k])
        while true
            y_best, j_best = -Inf, 0
            for j in ordering[1:k]
                if !has_edge(dag, j, i)
                    add_edge!(dag, j, i)
                    y′ = bayesian_score(n, r, dag, data, M, i)
                    if y′ > y_best
                        y_best, j_best = y′, j
                    end
                    rem_edge!(dag, j, i)
                end
            end
            if y_best > y
                y = y_best
                add_edge!(dag, j_best, i)
            else
                k==n-1 ? println("Bayesian score = $y_best") : print("")
                break
            end
        end
    end

    write_gph(dag, idx2names, outfile)
    return dag, idx2names
end

if length(ARGS) != 2
    error("usage: julia project1.jl <infile>.csv <outfile>.gph")
end

inputfilename = ARGS[1]
outputfilename = ARGS[2]

t_init = time()
dag, idx2names = compute(inputfilename, outputfilename)
t_f = time()
T = t_f-t_init
println("Done, time used = $T sec")

