clear,clc,close all
s = tf("s");
G = 0.1*exp(-6*s)/(s+0.1);
G = pade(G,4);
% bode(G)
% figure
% rlocus(G)
% figure
% step(4.33*G)

a = 0.6; b=0.5; g = 0.125;
Ku = 14.9;
omegau = 1.49;
Tu = 1/omegau;
Kp = a*Ku;
TI = b*Tu;
TD = g*Tu;
D = Kp*(1+1/TI/s+TD*s);
L = D*G;
figure
bode(L)
T = L/(1+L);
pStep(15,20,T,10)
%%
clear,clc,close all
z = tf("z",2);
G2 = 0.10136*(z+1)/z^4/(z-0.81594);
T2 = 1/z^4;
D2 = T2/G2/(1-T2);
pStep(10,35,T2,10)
az = @(m) m^4*(m-0.81594);
bz = @(m) 0.10136*(m+1);
D3 = az(z)/bz(1)/(z^5-bz(z)/bz(1));
T3 = G2*D3/(1+G2*D3);
pStep(10,35,T3,30)

function pStep(A,In,T,t)
    opt = stepDataOptions("StepAmplitude",A,"InputOffset",In);
    figure
    step(T,t,opt)
end