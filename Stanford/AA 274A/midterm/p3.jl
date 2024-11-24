using LinearAlgebra

th1 = 0:0.01:2*pi;
th2 = [pi+(pi/2-th-acos(4*sin(th)/2/sqrt(2))) for th in th1];

scatter(th1,th2)
