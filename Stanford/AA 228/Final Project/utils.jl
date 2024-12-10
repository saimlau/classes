using Distributions
using Base.Threads
using Combinatorics

Base.Dict{Symbol,V}(a::NamedTuple) where V =                         # Appendix G
    Dict{Symbol,V}(n=>v for (n,v) in zip(keys(a), values(a)))
Base.convert(::Type{Dict{Symbol,V}}, a::NamedTuple) where V =
    Dict{Symbol,V}(a)
Base.isequal(a::Dict{Symbol,<:Any}, nt::NamedTuple) =
    length(a) == length(nt) &&
    all(a[n] == v for (n,v) in zip(keys(nt), values(nt)))

struct SetCategorical{S}
    elements::Vector{S} # Set elements (could be repeated)
    distr::Categorical # Categorical distribution over set elements

    function SetCategorical(elements::AbstractVector{S}) where S
        weights = ones(length(elements))
        return new{S}(elements, Categorical(normalize(weights, 1)))
    end

    function SetCategorical(
        elements::AbstractVector{S},
        weights::AbstractVector{Float64}
        ) where S
        ℓ₁ = norm(weights,1)
        if ℓ₁ < 1e-6 || isinf(ℓ₁)
            return SetCategorical(elements)
        end
        distr = Categorical(normalize(weights, 1))
        return new{S}(elements, distr)
    end
end

Distributions.rand(D::SetCategorical) = D.elements[rand(D.distr)]
Distributions.rand(D::SetCategorical, n::Int) = D.elements[rand(D.distr, n)]
function Distributions.pdf(D::SetCategorical, x)
    sum(e == x ? w : 0.0 for (e,w) in zip(D.elements, D.distr.p))
end

####################################### 
#       Model-Specific stuff
#######################################
const DAYSLEFT = Dict(
    0 => -15,
    1 => -1,
    2 => 0,
    3 => 0,
    4 => 0,
    5 => 0,
    6 => 10
);
const MENTALSTATE = Dict(
    0 => -50,
    1 => 3,
    2 => 3,
    3 => 3,
);

function id2state(s)  # mapping from base 10 to base 7
    x = Int(floor(s/(16807)))
    tem = s%16807
    y = Int.(zeros(5))
    for i in eachindex(y)
        tem2 = 7^(5-i);
        y[i] = Int(floor(tem/tem2))
        tem %= tem2
    end
    return (x,y)
end

function state2id(x,y)
    out = x*16807;
    out += sum(y.*7.0.^(4:-1:0));
    return Int(out)
end

function id2obs(o)
    return id2state(o)
end

function obs2id(z,y)
    return state2id(z,y)
end

function possTransit(s,a)
    x,y = id2state(s)
    xs = Int.([])
    ys = []
    if a==6
        push!(xs,1,2,3)
        tem = copy(y)
        tem[y.!=6] .-= 1
        tem[tem.<0] .= 0
        doneIdx = findall(y.==6)
        push!(ys,tem)
        for i in combinations(doneIdx)
            tem2 = copy(tem)
            tem2[i] .= 5
            push!(ys,tem2)
        end
    else
        if x==0
            push!(xs,x)
            tem = copy(y)
            tem[a] = 6
            push!(ys,y,tem)
        else
            push!(xs,x,x-1)
            tem = copy(y)
            tem[a] = 6
            push!(ys,tem)
        end
    end
    return [state2id(xx,yy) for xx in xs for yy in ys]
end

function Transistion(s,a,s′)
    if !(s′ in possTransit(s,a))
        return 0.0
    end
    x,y = id2state(s)
    x′,y′ = id2state(s′)
    t = 1.0
    if a==6
        t *= 1/3
        for i in eachindex(y)
            if y[i]==6
                t *= y′[i]==5 ? 0.3 : 0.7
            end
        end
    else
        if x==0
            t *= y′[a]==6 ? 0.2 : 0.8
        else
            t *= x′==x ? 0.2 : 0.8
        end
    end
    return t
end

function possPreTransit(a,s′)
    x′,y′ = id2state(s′)
    xs = Int.([])
    ys = []
    if a==6
        push!(xs,0,1,2,3)
        tem = copy(y′)
        tem[y′.<=4] .+= 1
        push!(ys,tem)
        for i in combinations(findall(y′.==0))
            tem3 = copy(tem)
            tem3[i] .= 0
            push!(ys,tem3)
        end
        doneIdx = findall(y′.==5)
        tem4 = copy(ys)
        for tem5 in tem4
            for i in combinations(doneIdx)
                tem2 = copy(tem5)
                tem2[i] .= 6
                push!(ys,tem2)
            end
        end
    else
        if x′==0
            push!(xs,0,1)
            push!(ys,y′)
            if y′[a]==6
                for j in 0:6
                    tem = copy(y′)
                    tem[a] = j
                    push!(ys,tem)
                end
            end
        else
            push!(xs,x′,x′+1)
            for j in 0:6
                tem = copy(y′)
                tem[a] = j
                push!(ys,tem)
            end
        end
    end
    return [state2id(xx,yy) for xx in xs for yy in ys]
end

function Reward(s,a)
    x,y = id2state(s)
    r = MENTALSTATE[x]
    r += sum([DAYSLEFT[yi] for yi in y])
    r += a==6 ? 0 : -3
    return r
end

function possObs(a,s′)
    x′,y′ = id2state(s′)
    zs = []
    if x′==0
        if a==6
            push!(zs,0,1)
        else
            push!(zs,0)
        end
    else
        push!(zs,x′)
        if x′ in [1,2]
            push!(zs,x′-1,x′+1)
        else
            push!(zs,2)
        end
    end
    return [obs2id(zz,y′) for zz in zs]
end

function ObsModel(a, s′, o)
    if !(o in possObs(a,s′))
        return 0.0
    end
    x′, ~ = id2state(s′)
    z, ~ = id2obs(o)
    p = 1.0
    if x′==0
        if a==6
            p *= z==0 ? 0.7 : 0.3
        end
    else
        if x′==3
            p *= z==2 ? 0.4 : 0.6
        else
            p *= z==x′ ? 0.6 : 0.2
        end
    end
    return p
end

dis2 = SetCategorical([6,5],[0.7,0.3]);
dis3 = SetCategorical([-1,0],[0.8,0.2]);
function TRO(s,a)                        # x in {0,1,2,3}
    x,y = id2state(s)
    r = Reward(s,a)
    if a ==6 || (x==0 && rand(SetCategorical([0,1],[0.5,0.5]))==1)
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
    o = observe(a,x′,y′)
    # if !(state2id(x′,y′) in possTransit(s,a))
    #     println(x′)
    #     println(y′)
    #     error("Getting impossible state!!!!!")
    # end
    return (state2id(x′,y′), r, o)
end

dis4 = SetCategorical([0,1],[0.7,0.3]);
dis5 = SetCategorical([2,3],[0.4,0.6]);
dis6 = SetCategorical([-1,0,1],[0.2,0.6,0.2]);
function observe(a,x,y)
    z = -1
    if x==0
        if a==6
            z = rand(dis4)
        else
            z = 0
        end
    else
        if x==3
            z = rand(dis5)
        else
            z = x+rand(dis6)
        end
    end
    return obs2id(z,y)
end

function TR(s,a)                        # x in {0,1,2,3}
    x,y = id2state(s)
    r = Reward(s,a)
    if a ==6 || (x==0 && rand(SetCategorical([0,1],[0.7,0.3]))==1)
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
#######################################

struct POMDP
    γ   # discount factor
    S   # state space
    A   # action space
    𝒪   # observation space
    T   # transition function
    R   # reward function
    O   # observation function
    TRO # sample transition, reward, and observation
end

struct AlphaVectorPolicy    # Algorithm 20.4
    P # POMDP problem
    Γ # alpha vectors
    a # actions associated with alpha vectors
end
function utility(π::AlphaVectorPolicy, b)
    return maximum(α⋅b for α in π.Γ)
end
function (π::AlphaVectorPolicy)(b)
    i = argmax([α⋅b for α in π.Γ])
    return π.a[i]
end

function alphavector_iteration(P::POMDP, M, Γ)  # Algorithm 21.1
    for k in 1:M.k_max
        Γ = update(P, M, Γ)
    end
    return Γ
end

struct QMDP   # Algorithm 21.2
    k_max # maximum number of iterations
end

function update(P::POMDP, M::QMDP, Γ)
    S, A, R, T, γ = P.S, P.A, P.R, P.T, P.γ
    Γ′ = [[R(s,a) + γ*sum(T(s,a,s′)*maximum(α′[j] for α′ in Γ)
        for (j,s′) in enumerate(possTransit(s,a))) for s in S] for a in A]
    return Γ′
end

function solve(M::QMDP, P::POMDP)
    Γ = [zeros(length(P.S)) for a in P.A]    # Initiates 1 all zeros α vector for each action
    Γ = alphavector_iteration(P, M, Γ)
    return AlphaVectorPolicy(P, Γ, P.A)
end

function baws_lowerbound(P::POMDP)        # Algorithm 21.4, Eqn. 21.5
    S, A, R, γ = P.S, P.A, P.R, P.γ
    r = maximum(minimum(R(s, a) for s in S) for a in A) / (1-γ)
    α = fill(r, length(S))
    return α
end

struct MDP           # Algorithm 7.1
    γ # discount factor
    S # state space
    A # action space
    T # transition function
    R # reward function
    TR # sample transition and reward
end

function plotResults(name,T,rr,xx,yy)
    plot(1:T,rr, label="", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=1, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottomright, background_color = :transparent)
    xlabel!("time (days)")
    ylabel!("Reward")
    # title!(name)
    savefig("rewardOver$(T)Days($name).png")

    acR = accumulate(+,rr);
    plot(1:T,acR, label="", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=1, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottomright, background_color = :transparent)
    xlabel!("time (days)")
    ylabel!("accumulated utility")
    # title!(name)
    savefig("accRewardOver$(T)Days($name).png")

    plot(0:T,xx, label="", framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=1, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottomright, background_color = :transparent)
    xlabel!("time (days)")
    ylabel!("Mental State")
    # title!(name)
    savefig("MStateOver$(T)Days($name).png")

    plot(framestyle = :box, xguidefontsize=12, yguidefontsize=12,legendfontsize=12, ytickfontsize = 12, xtickfontsize = 12,
    linewidth=1, background_color_legend = nothing, foreground_color_legend = nothing, legend=:bottomright, background_color = :transparent)
    for i in 1:5
        plot!(0:T,yy[:,i], label="Task $i")
    end
    xlabel!("time (days)")
    ylabel!("Days Until Due")
    # title!(name)
    savefig("TStateOver$(T)Days($name).png")
end