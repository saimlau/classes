clear,clc,close all
Ds = RR_tf([-1],[0 -10],1);
%Ds = RR_tf([],[-2 -1],3);
% Dz = SML_C2D_matched(Ds,1,0,true)
Dz = SML_C2D_matched(Ds,1,0.0001i,true)

Dzzz = c2d(zpk([-1],[0 -10],1),1,'matched')

syms z1 p1 omega
Ds2 = RR_tf([-z1],[0 -p1],1);
Dz2 = SML_C2D_matched(Ds2,1,omega,true);