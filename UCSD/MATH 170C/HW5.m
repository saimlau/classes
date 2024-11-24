clear,clc,close all
f = @(t,y,yp) exp(t)+y*cos(t)-(t+1)*yp;
alpha = 1; beta = 3; a = 0; b = 1;
hs=[0.1,0.05]; mkr = ["b-","g--"];
for i=1:2
    h = hs(i);
    [x1,t] = BVP_shooting(f,a,b,alpha,beta, h, 0.0001, 500);
    plot(t,x1, mkr(i), LineWidth=2), hold on
end

mkr = ["m:", "r*"];
for i=1:2
    h = hs(i);
    [x1,t] = BVP_finitediff(@(t) exp(t),@(t) cos(t),@(t) -(t+1),a,b,alpha,beta, h);
    plot(t,x1, mkr(i), LineWidth=2, MarkerSize=6), hold on
end
grid on
legend("h=0.1, Shooting", "h=0.05, Shooting", "h=0.1, Finitediff", "h=0.05, Finitediff", ...
    Location="southeast")
set(gca, "LineWidth", 2,"FontSize", 12)