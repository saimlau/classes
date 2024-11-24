# problem 2
include("readclassjson.jl");
using LinearAlgebra
using Plots

data = readclassjson("HW2/color_perception_data.json");
test_light = data["test_light"];
wavelen = data["wavelength"];
tungsten = data["tungsten"];
sunlight = data["sunlight"];
l = data["L_coefficients"];
m = data["M_coefficients"];
s = data["S_coefficients"];
R = data["R_phosphor"];
G = data["G_phosphor"];
B = data["B_phosphor"];
n = 20;

# Part c
T1 = Transpose([l m s]);
Goal = T1*test_light;
T2 = T1*[R G B];
a_sol = inv(T2)*Goal;
println("a_sol = $a_sol")

# Part d
r1 = zeros(n);
r2 = zeros(n);
r1[1:3] = [0.3 0.5 0.6];
p1 = r1.*tungsten;
cone_1 = T1*p1;
p2 = T1[1:3,2:4]\cone_1;
r2[2:4] = p2./tungsten[2:4];
println("r1 = $r1")
println("r2 = $r2")

T1*(r1.*tungsten)
T1*(r2.*tungsten)
T1*(r1.*sunlight)
T1*(r2.*sunlight)

