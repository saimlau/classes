using LinearAlgebra
print("\033c")
### Fully Observable
#Q3
b = [0,-1,0];
tem = 0.9*[0.7 0.3 0;0 1 0;0 0 1];
U3 = (I-tem)\b
println("U3=$U3")

#Q4
b = [1.1,0,0];
tem = 0.9*[1 0 0;0.5 0.5 0;0 0 1];
U4 = (I-tem)\b
println("U4=$U4")

#Q5
A = 1:3;
S = 1:3;
γ = 0.9;
R = zeros(length(S),length(A));
R[1,2] = 1.1;
R[1,3] = 10;
R[2,1] = -1;
T = zeros(length(S),length(A),length(S));
T[1,1,1] = 0.7;
T[1,1,2] = 0.3;
T[1,2,1] = 1;
T[2,1,2] = 1;
T[2,2,1] = 0.5;
T[2,2,2] = 0.5;
T[:,3,3] .= 1;
T[3,:,3] .= 1;
U = zeros(length(S));
i = 0;
while true
    U_prev = copy(U)
    for s in S
        U[s] = maximum([R[s,a] + γ*dot(T[s,a,:],U_prev) for a in A])
    end
    i +=1
    if norm(U_prev-U)<1e-12
        break
    end
end
println("Converged in $i iterations,")
println("U*=$U")

#Q6
sol = argmax([R[2,a]+γ*dot(T[2,a,:],U) for a in A]);
println("π*(sₙ)=$sol")

#Q7
sol = argmax([R[1,a]+γ*dot(T[1,a,:],U) for a in A]);
println("π*(sᵢ)=$sol")

### Partially Observable
#Q8
println("MDP->POMDP")
𝒪 = 1:3;
function update(b::Vector{Float64}, S, T, O, a, o)
    b′ = similar(b)
    for (i′, s′) in enumerate(S)
        po = O[s′, o]
        b′[i′] = po * sum(T[s, a, s′] * b[i] for (i, s) in enumerate(S))
    end
        if sum(b′) ≈ 0.0
        fill!(b′, 1)
    end
    return normalize!(b′, 1)
end
O = zeros(length(S),length(𝒪));
O[1,1] = 0.7;
O[1,2] = 0.05;
O[1,3] = 0.25;
O[2,1] = 0.05;
O[2,2] = 0.25;
O[2,3] = 0.7;
O[3,3] = 1;
b = [0.5,0.5,0];
b = update(b,S,T,O,1,3);
println("b8′=$b")

#Q9
b = [1.0,0,0];
for i in 1:3
    b = update(b,S,T,O,1,1);
end
println("b9′=$b")

#Q15
r = maximum(minimum(R[s, a] for s in S) for a in A) / (1-γ);
α = fill(r, length(S));
println("αBAWS=$α")

#Q16
function BLB(α, S, R, T, γ, a, iters)
    α_new = similar(α)
    for i in 1:iters
        for s in S
            α_new[s] = R[s,a]+γ*dot(T[s,a,:],α)
        end
        α = copy(α_new)
    end
    return α_new
end
α = [-10.0,-10.0,-10.0];
αT = BLB(α,S,R,T,γ,1,2)
αA = BLB(α,S,R,T,γ,2,2)
αE = BLB(α,S,R,T,γ,3,2)
println("αT=$αT")
println("αA=$αA")
println("αE=$αE")

