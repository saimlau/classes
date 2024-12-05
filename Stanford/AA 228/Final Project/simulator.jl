include("utils.jl");

function updateb(b::Vector{Float64}, P, a, o)         # 
    S, T, O = P.S, P.T, P.O
    b′ = similar(b)
    for (i′, s′) in enumerate(S)
        po = O(a, s′, o)
        if po==0.0
            continue
        end
        b′[i′] = po * sum(T(s, a, s′) * b[i] for (i, s) in enumerate(possPreTransit(a,s′)))
    end
    if sum(b′) ≈ 0.0
        fill!(b′, 1)
    end
    if any(isnan.(b))
        error("Nan in b !!!")
    end
    return normalize!(b′, 1)
end;