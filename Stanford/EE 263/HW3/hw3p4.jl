# Problem 4
include("readclassjson.jl");
using LinearAlgebra
print("\033c")
data = readclassjson("HW3/lti_memory_data.json");
T = data["T"];
u = data["u"];
y = data["y"];

M_max = Int(ceil(T/2))

function test_good(h)
    M_test = length(h)
    U_test = zeros(M_test,M_test)
    Y_test = zeros(M_test)
    T_test = T-5
    for i =1:M_test
        U_test[M_test+1-i,:] = u[Vector((T_test-i):-1:(T_test-i+1-M_test))]
        Y_test[M_test+1-i] = y[T_test-i+1]
    end
    return norm(U_test*h-Y_test)<0.05
end

for M=2:M_max
    global U = zeros(M,M)
    global Y = zeros(M)

    for i =1:M
        U[M+1-i,:] = u[Vector((T-i):-1:(T-i+1-M))]
        Y[M+1-i] = y[T-i+1]
    end
    try
        global h = inv(U)*Y
        if test_good(h)
            println("Soln found")
            global M_sol = M
            println("M = $M")
            break
        end
    catch
        println("singular, $M")
    end
end
