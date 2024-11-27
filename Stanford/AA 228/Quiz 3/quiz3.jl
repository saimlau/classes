using LinearAlgebra
include("utils.jl")


function update(b::Vector{Float64}, P, a, o)
    S, T, O = P.S, P.T, P.O
    b′ = similar(b)
    for (i′, s′) in enumerate(S)
        po = O(a, s′, o)
        b′[i′] = po * sum(T(s, a, s′) * b[i] for (i, s) in enumerate(S))
    end
        if sum(b′) ≈ 0.0
        fill!(b′, 1)
    end
    return normalize!(b′, 1)
end