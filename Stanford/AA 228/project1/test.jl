using LinearAlgebra

print("\033c")
function update_M(M, i)
    M[i] = ones(3,5)
end

function gogo(n)
    M = [zeros(i,i+1) for i in 1:n]
    println(M)

    for i in 2:n-1
        update_M(M,i)
    end
    for i in 1:n
        println(M[i])
    end
end

gogo(6)