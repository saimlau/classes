%2b.  Is the controller determined in 1a proper?  If not, change your target f(s) by adding 
% k additional poles to T(s) at, say, s=-20, for a sufficiently large k such that the answer 
% returned by RR_Diophantine is proper, compute the modified D(s), and check that it worked.  
% (How big does k need to be?). Discuss.
clear,clc,close all
p = ones(1,6)*-20;
a = RR_poly([-1 1 -3 3 -6 6],1); b = RR_poly([-2 2 -5 5],1); f = RR_poly([-1,-1,-3,-3,-6,-6, p],1); G = RR_tf(b,a);
[x2,y2] = RR_diophantine(a,b,f); D = RR_tf(y2,x2)
test = trim(a*x2+b*y2), residual = norm(trim(a*x2+b*y2)-f)