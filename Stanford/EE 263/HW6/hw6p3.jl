# Problem 3
using LinearAlgebra
using Plots
using LaTeXStrings
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW6/backwards_compatible_transceiver_data.json");
Anew = data["Anew"]
m = data["m"]
Aold = data["Aold"]
x = data["x"]
n = data["n"]

# Part b)
function transformA(A)
    AA = vcat([hcat([transpose(A).*j for j in I(n)[i,:]]...) for i in 1:n]...)
    return AA
end
Atnew = transformA(Anew);
Atold = transformA(Aold);
vI = vec(I(n));
F = [2*transpose(Atold)*Atold transpose(Atnew);Atnew zeros(n^2,n^2)];
dd = [2*transpose(Atold)*vI;vI];
vB = (inv(F)*dd)[1:n*m];
B = transpose(reshape(vB,(m,n)));
J1 = norm(B*Aold-I)^2
println("Optimal J = $J1")
J2 = norm(B*Anew-I)^2

B2 = inv(transpose(Anew)*Anew)*transpose(Anew);
J3 = norm(B2*Aold-I)^2
println("J from A† = $J3")
J4 = norm(B2*Anew-I)^2

# Part c)
yold = Aold*x;
x_est = B*yold;
x_est = [xi>0.5 ? 1 : 0 for xi in x_est];
x_est2 = B2*yold;
x_est2 = [xi>0.5 ? 1 : 0 for xi in x_est2];
bit_err = sum(x_est.!=x)/n;
bit_err2 = sum(x_est2.!=x)/n;
println("Bit error rate = $bit_err")
println("Bit error rate from A† = $bit_err2")

