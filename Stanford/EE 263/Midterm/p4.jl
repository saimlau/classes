# Problem 4
using LinearAlgebra
using Plots
using LaTeXStrings
include("readclassjson.jl");
print("\033c")
data = readclassjson("Midterm/power_system.json");
V_sm = data["V_small"]
I_sm = data["I_small"]
Y_bus_sm = data["Y_bus_small"]
V_l1 = data["V_large_1"]
V_l2 = data["V_large_2"]
I_l1 = data["I_large_1"]
I_l2 = data["I_large_2"]

function form_A(V)
    n, M = size(V)
    r = Int((n^2-n)/2)
    A_f = zeros((0,r));
    for m in 1:M
        A = zeros((n,r))
        for k in 1:n
            q = 1
            for i in 1:n-1
                for j in i+1:n
                    if i==k
                        A[k,q] = V[k,m]-V[j,m]
                    elseif j==k
                        A[k,q] = V[k,m]-V[i,m]
                    end
                    q+=1;
                end
            end
        end
        A_f = vcat(A_f,A)
    end
    return A_f
end;
function form_I(I_da)
    n,M = size(I_da)
    I_f = I_da[:,1]
    for m in 2:M
        I_f = vcat(I_f,I_da[:,m])
    end
    return I_f
end;

# Part c
function solve(V_dat, I_dat, ls)
    n, M = size(V_dat)
    r = Int((n^2-n)/2)
    A = form_A(V_dat);
    I_stk = form_I(I_dat);
    if ls
        y = inv(transpose(A)*A)*transpose(A)*I_stk;
        println("rank(A) = $(rank(A))")
    else
        # y = transpose(A)*inv(A*transpose(A))*I_stk;
        println("rank(A) = $(rank(A))")
        println("size(A) = $(size(A))")
        a = nullspace(A)[:,1]
        b = nullspace(A)[:,2]
        c = nullspace(A)[:,3]
        y = inv(transpose(A)*A)*transpose(A)*I_stk-a.*900-b.*600+c.*300;
    end
    Y = zeros((n,n));
    k = 1;
    for i in 1:n
        for j in i+1:n
            Y[i,j] = -y[k]
            Y[j,i] = -y[k]
            k += 1
        end
    end
    Y[diagind(Y)] = -sum(Y,dims=2);
    Y_cut = copy(Y);
    Y_cut[abs.(Y).<20] .= 0.0;
    k = r;
    for i in n:-1:1
        for j in n:-1:i+1
            if Y_cut[i,j] == 0.0
                A = A[:,1:size(A)[2] .!= k]
            end
            k -= 1;
        end
    end
    if ls
        y_sp = inv(transpose(A)*A)*transpose(A)*I_stk;
    else
        # y_sp = transpose(A)*inv(A*transpose(A))*I_stk;
        y_sp = A\I_stk;
    end
    Y_sp = zeros((n,n));
    k = 1;
    for i in 1:n
        for j in i+1:n
            if Y_cut[i,j] != 0
                Y_sp[i,j] = -y_sp[k]
                Y_sp[j,i] = -y_sp[k]
                k += 1
            end
        end
    end
    Y_sp[diagind(Y_sp)] = -sum(Y_sp,dims=2);
    return Y, Y_sp
end

# Part c
Y, Y_sp = solve(V_sm,I_sm,true)
print("Y^(BLUE) = ")
display(Y)
println("||Y^(BLUE)-Y_bus|| = $(norm(Y-Y_bus_sm))")
print("\nY^(SPARSE-BLUE) = ")
display(Y_sp)
print("||Y^(SPARSE-BLUE)-Y_bus|| = $(norm(Y_sp-Y_bus_sm))")
println(" ")

# Part d
Y_l1, Y_sp_l1 = solve(V_l1,I_l1,true);
n = size(Y_l1)[1];
edge_cnt = (n^2-n)/2-length(findall(Y_sp_l1.==0))/2;
println("\nEdges count (Large1) = $edge_cnt")
I_est_l1 = Y_sp_l1*V_l1;
println("I error (Large 1) = $(norm(I_l1-I_est_l1))")

# Part e
Y_l2, Y_sp_l2 = solve(V_l2,I_l2,false);
print("Y^(BLUE, adjusted) (Large 2) = ")
display(Y_l2)