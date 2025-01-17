clear,clc,close all
ths = 2*pi/3:-0.01:pi/3;
L = 1;
m = 1;
g = 9.81;
thd0 = -1;
x = @(th) L.*cos(th);
y = @(th) L.*sin(th);
figure
plot(ths,y(ths),LineWidth=2)
xlabel("$\theta$",Interpreter="latex")
ylabel("$y\;(m)$",Interpreter="latex")
set(gca,"LineWidth",2,"FontSize",15,"XDir","reverse")
xticks(pi/3:pi/3:2*pi/3)
xticklabels(["π/3","2π/3"])
ylim([0,1.2])
figure
plot(ths,-m*L.*sin(ths)+m*g,LineWidth=2)
xlabel("$\theta$",Interpreter="latex")
ylabel("$N_y\;(N)$",Interpreter="latex")
set(gca,"LineWidth",2,"FontSize",15,"XDir","reverse")
xticks(pi/3:pi/3:2*pi/3)
xticklabels(["π/3","2π/3"])
figure
hold on

ts = 0:0.001:(-1+sqrt(1+2*pi/3));
thdd = -1;
thd = -1-ts;
th = 2*pi/3-ts-ts.^2./2;
Ny = m*L*(thdd.*cos(th)-thd.^2.*sin(th))+m*g;

plot(th,Ny,LineWidth=2)
xlabel("$\theta$",Interpreter="latex")
ylabel("$N_y\;(N)$",Interpreter="latex")
set(gca,"LineWidth",2,"FontSize",15,"XDir","reverse","Box","on")
xticks(pi/3:pi/3:2*pi/3)
xticklabels(["π/3","2π/3"])