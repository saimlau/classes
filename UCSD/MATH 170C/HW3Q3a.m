clear,clc,close all
x0 = 0; a = 0; b = 5; h = 0.01;
lamdas = [5,-5,-10];
maks = ["b.","r*","g--"];
for i=1:3
    lamda = lamdas(i);
    f = @(t,x) lamda.*x+cos(t)-lamda.*sin(t);
    [x,t] = RK4(f,[x0],a,b,h);
    plot(t,x,maks(i),MarkerSize=6,LineWidth=3), hold on
end
plot(t,sin(t),'k-',LineWidth=2)
legend('\lambda=5','\lambda=-5','\lambda=-10','Exact',Location='northwest')
set(gca,"LineWidth",2,"FontSize",12)
xlabel('t'), ylabel('x')
grid on