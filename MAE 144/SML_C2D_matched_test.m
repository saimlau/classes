% 2 tests for SML_C2D_matched, Ds(s) = (s+z1)/(s(s+p1)), 
% where z1=1 and p1=10 for the first test, and symbolic for the second.
clear,close all
disp("Test 1: Comparing with MATLAB provided c2d")
Ds = RR_tf([-1],[0 -10],1)
h = 1; omega = 0.000001i;
Dz_SML_C2D = SML_C2D_matched(Ds,h,omega,true)
Dz_Matlab = c2d(zpk([-1],[0 -10],1),1,'matched')
Dss = tf(Ds.num.poly,Ds.den.poly);
Dzz = tf(Dz_SML_C2D.num.poly,Dz_SML_C2D.den.poly,h); 
bode(Dss,Dz_Matlab,Dzz),legend

pause
disp("Test 2: Symbolic calculation, compare with pen and paper calculation.")
syms z1 p1 omega
Ds2 = RR_tf([-z1],[0 -p1],1)
Dz2 = SML_C2D_matched(Ds2,1,omega,true)