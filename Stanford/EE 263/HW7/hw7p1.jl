# Problem 1
using LinearAlgebra
print("\033c")

Am = [2 1.2 -1;0.4 2 -0.5; -0.5 0.9 1];
v1 = [0.7,0,0.7];
v2 = [0.3,0.6,0.7];
v3 = [0.6,0.6,0.3];
T = [v1 v2 v3];
U = inv(T);
n = size(Am)[1];

B = zeros(n^2,n);
for i in 1:n
    for j in 1:n
        B[(j-1)*n+i,:] = [T[i,k]*U[k,j] for k in 1:n]
    end
end

lT = inv(transpose(B)*B)*transpose(B)*Am[:];
als = B*lT;
Â = reshape(als,n,n);
J = norm(Am[:]-als)^2/n^2;
println("Â = ")
display(Â)
println("J for Â: $J")
println("Verification:")
println("Âv₁./v₁ = $(Â*v1./v1)")
println("Âv₂./v₂ = $(Â*v2./v2)")
println("Âv₃./v₃ = $(Â*v3./v3)")
println("λ=$lT")
