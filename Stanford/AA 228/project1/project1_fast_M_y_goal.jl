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

Random.seed!(87577658735718327237845738)
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
    draw(PDF(fname*"3.pdf", 16cm, 16cm), p)
end

function prior(n::Int, r::Vector{Int64}, G::DiGraph)
    q = [prod([r[j] for j in inneighbors(G,i)]) for i in 1:n]
    return [ones(q[i], r[i]) for i in 1:n]
end

function init_M(n, r, G, D)
    q = [prod([r[j] for j in inneighbors(G,i)]) for i in 1:n]
    M = [zeros(q[i], r[i]) for i in 1:n]
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

function sub2ind(siz, x)
    k = vcat(1, cumprod(siz[1:end-1]))
    return dot(k, x .- 1) + 1
end

function update_M(r, G, D::DataFrame, M, i_affected)
    q_i = prod([r[j] for j in inneighbors(G,i_affected)])
    M[i_affected] = zeros(q_i, r[i_affected])

    for o in eachrow(D)
        k = o[i_affected]
        parents = inneighbors(G,i_affected)
        j = 1
        if !isempty(parents)
            j = sub2ind(r[parents], Vector(o[parents]))
        end
        M[i_affected][j,k] += 1.0
    end
end

function bayesian_score_component(M, A)
    p = sum(loggamma.(A + M))
    p -= sum(loggamma.(A))
    p += sum(loggamma.(sum(A,dims=2)))
    p -= sum(loggamma.(sum(A,dims=2) + sum(M,dims=2)))
    return p
end

function bayesian_score(n::Int, r::Vector{Int64}, G::DiGraph, D::DataFrame, M, i_affected)
    update_M(r, G, D, M, i_affected)
    alpha = prior(n, r, G)
    return sum(bayesian_score_component(M[i], alpha[i]) for i in 1:n)
end

function rand_graph_neighbor(G)
    n = nv(G)
    i = rand(1:n)
    j = mod1(i + rand(2:n)-1, n)
    G′ = copy(G)
    p = rand()
    if p>=1.0
        ran_iter = rand(1:n)
        i_affected = []
        for k in 1:ran_iter
            i = rand(1:n)
            j = mod1(i + rand(2:n)-1, n)
            has_edge(G, i, j) ? rem_edge!(G′, i, j) : add_edge!(G′, i, j)
            push!(i_affected, i, j)
        end
        return G′, unique(i_affected)
    elseif p>=0.0 && has_edge(G, i, j)
        rem_edge!(G′, i, j)
        add_edge!(G′, j, i)
    else
        has_edge(G, i, j) ? rem_edge!(G′, i, j) : add_edge!(G′, i, j)
    end
    return G′, [i j]
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
    
    ordering_k2_best = 1:n
    M = init_M(n,r,dag,data)
    M_initial = [copy(m) for m in M]
    y = bayesian_score(n, r, dag, data, M, 1)
    y_best_k2 = y
    Max_k2_itr = 1
    
    for l in 1:Max_k2_itr
        dag_tem = DiGraph()
        add_vertices!(dag_tem, n)
        M_tem = [copy(m) for m in M_initial]
        # ordering = shuffle(shuffle(1:n))
        ordering = [36, 26, 7, 42, 4, 29, 28, 27, 37, 50, 14, 32, 18, 48, 10, 35, 12, 23, 15, 20, 17, 39, 44, 47, 9, 46, 31, 13, 24, 43, 45, 16, 38, 8, 49, 11, 21, 40, 25, 41, 3, 6, 5, 1, 30, 22, 33, 34, 19, 2]
        for (k,i) in enumerate(ordering[2:end])
            y = bayesian_score(n, r, dag_tem, data, M_tem, ordering[k])
            while true
                y_best, j_best = -Inf, 0
                for j in ordering[1:k]
                    if !has_edge(dag_tem, j, i)
                        add_edge!(dag_tem, j, i)
                        y′ = bayesian_score(n, r, dag_tem, data, M_tem, i)
                        if y′ > y_best
                            y_best, j_best = y′, j
                        end
                        rem_edge!(dag_tem, j, i)
                    end
                end
                if y_best > y
                    y = y_best
                    add_edge!(dag_tem, j_best, i)
                else
                    if k==n-1 && y > y_best_k2
                        y_best_k2 = y
                        dag = copy(dag_tem)
                        ordering_k2_best = ordering
                        println(ordering)
                        println("Bayesian score = $y_best_k2")
                    end
                    break
                end
            end
        end
    end
    y = y_best_k2
    ordering = ordering_k2_best
    M = init_M(n,r,dag,data)

    # y = bayesian_score(n, r, dag, data, M, ordering[end])
    Max_itr = 5000
    if outfile=="medium.gph"
        y_goal = -96285.22362037576
    elseif outfile=="small.gph"
        y_goal = -3794.855597709796
    else
        # y_goal = -398870.6799705375
        y_goal = -399999
    end
    for k in 1:Max_itr
    # while y<=y_goal
        G′, i_changed = rand_graph_neighbor(dag)
        [update_M(r, G′, data, M, i_changed[i]) for i in 1:length(i_changed)-1]
        y′ = is_cyclic(G′) ? -Inf : bayesian_score(n, r, G′, data, M, i_changed[2])
        if y′ > y
            y, dag = y′[end], G′
        else
            for i in i_changed
                update_M(r, dag, data, M, i)
            end
        end
    end
    println("Final Bayesian score = $y")

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
