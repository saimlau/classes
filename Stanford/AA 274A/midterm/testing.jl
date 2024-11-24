using LinearAlgebra
using Symbolics  
using Rotations

@variables x

x = 5

# th = 5*x
th = deg2rad(5*x)
K = [100 0 640;0 100 480;0 0 1];
Rt = [cos(th) -sin(th) 0 x;sin(th) cos(th) 0 x;0 0 0 0];

P = Rt*[0.5;0.5;0;1]
P_tar = [5;5]
norm(P[1:2]-P_tar)

p = [400;340;1];
uv = inv(K)*p

K*Rt

# 2 iv
R1 = hcat(vcat(RotXYZ(0,deg2rad(45),0),[0 0 0]),[0;0;0;1])
R2 = hcat(vcat(RotXYZ(deg2rad(60),0,0),[0 0 0]),[0;0;0;1])
R = R2*R1
x = inv(R)*[1;1;1;1]