# Problem5
using LinearAlgebra
print("\033c")
# Part a
p=6;
q=11;
S = [[1 2 3],[1 2 4],[1 2 6],[1 3 5],[1 4 5],[2 3 6],[2 4 6],[3 4 5],[3 5 6],[4 5 6],[1 2 3 4 5 6]]
A = zeros((q,p))
v = [-2;14;6;4;20;-5;11;9;1;17;15]
for i in 1:q
    for j in S[i]
        A[i,j] = 1.0
    end
end
println("null(A) = ")
println(nullspace(A))
println("rank(A) = ")
println(rank(A))
u = A\v;


u = round.(u,sigdigits=5)
println("u = \n $u")
println("Au1 = ")
println(round.(A*u))

if !isempty(nullspace(A))
    u2 = u+nullspace(A)
    println("u2 = \n $u2")
    println("Au2 = ")
    println(round.(A*u2))
end