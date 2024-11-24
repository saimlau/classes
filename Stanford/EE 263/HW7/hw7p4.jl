# Problem 4
using LinearAlgebra
using LaTeXStrings
include("readclassjson.jl");

print("\033c")
data = readclassjson("HW7/two_firms_data.json");
B = data["B"]
c = data["c"]
a = data["a"]

q1 = inv(B)*(a-c)./2;
q2 = inv(B)*(a-B*q1-c)./2;
p = a-B*(q1+q2);
π1 = transpose(p-c)*q1;
π2 = transpose(p-c)*q2;
println("q⁽¹⁾=$(round.(q1,sigdigits=4))")
println("q⁽²⁾=$(round.(q2,sigdigits=4))")
println("p=$(round.(p,sigdigits=4))")
println("π₁=$(round.(π1,sigdigits=6))")
println("π₂=$(round.(π2,sigdigits=6))")