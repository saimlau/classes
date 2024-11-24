# problem 6
include("readclassjson.jl")
using LinearAlgebra
using Plots

data = readclassjson("HW2/robot_coin_collector.json");
n = data["n"];
x = data["x"];

# A = LowerTriangular(4*ones(n,n));
# A[diagind(A)] = 2*ones(n);
A = [1 0 0 0 0 0; 3 1 0 0 0 0; 5 3 1 0 0 0; 7 5 3 1 0 0; 9 7 5 3 1 0; 11 9 7 5 3 1].*2.0;
B = [0 0 0 0 0 0; 1 0 0 0 0 0; 2 1 0 0 0 0; 3 2 1 0 0 0; 4 3 2 1 0 0; 5 4 3 2 1 0].*4.0;
B[diagind(B)] = 0.5*ones(n);
x_even = x[2:2:2*n];
f = A\x_even;

x_odd = x[1:2:2*n];
# B = LowerTriangular(2*ones(n,n));
# B[diagind(B)] = 1/2*ones(n);
x_odd_til = B*f;

x_odd_til.==x_odd

C = vcat([A;B]);
x_new = vcat(x_even, x_odd);

for i=1:2*n
    C_til = vcat([C[1:i-1,:];C[i+1:2*n,:]])
    x_new_til = vcat(x_new[1:i-1], x_new[i+1:2*n])
    f_til = C_til\x_new_til
    if norm(C_til*f_til-x_new_til) < 2
        println("Success")
        println("f = $f_til")
        global f_sol = f_til
        i<n ? coin = 2*i : coin = (i-n)*2-1
        println("Coin $coin not collected.")
        break
    end
end

A_final = vcat([B[1:1,:];A[1:1,:]])
for i=2:n
    tem = vcat([B[i:i,:];A[i:i,:]]);
    A_final = vcat([A_final; tem])
end
x_path = A_final*f;
pushfirst!(x_path, 0)

r = scatter(x_path, 0:2*n, label="Robot Traj.", reuse=true, xlabel="x", ylabel="y")
c = scatter!(x, 1:2*n, label="Coins",markershape=:star5)
savefig("HW2/hw2p6.png")
