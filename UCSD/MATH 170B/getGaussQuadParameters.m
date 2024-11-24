clear,clc,close all
syms c0 c1 c2 c3 x0 x1 x2 x3
ff = @(n) c0*x0^n+c1*x1^n+c2*x2^n+c3*x3^n;
f0 = 2==ff(0);
f1 = 0==ff(1);
f2 = 2/3==ff(2);
f3 = 0==ff(3);
f4 = 2/5==ff(4);
f5 = 0==ff(5);
f6 = 2/7==ff(6);
f7 = 0==ff(7);
f8 = x0<x1;
f9 = x1<x2;
f10 = x2<x3;

sol = solve(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,c0,c1,c2,c3,x0,x1,x2,x3);
%%
c0 = simplify(sol.c0);
c1 = simplify(sol.c1);
c2 = simplify(sol.c2);
c3 = simplify(sol.c3);
x0 = simplify(sol.x0);
x1 = simplify(sol.x1);
x2 = simplify(sol.x2);
x3 = simplify(sol.x3);
save('4ptGaussQuadParameters.mat')