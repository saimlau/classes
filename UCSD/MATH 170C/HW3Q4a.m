clear,clc,close all
x0 = 1; a = 0; b = 1; h = 0.25; TOL = 0.001; MaxIters = 50;
f = @(t,x) -2.*t.*x.^2;
[x,t] = AM4(f,x0,a,b,h,TOL,MaxIters);
x_exac = 1./(1+t.^2);
plot(t,x_exac,'g-',LineWidth=2), hold on
plot(t,x,'b--',LineWidth=2)
legend('Exact','AM4')
set(gca,"LineWidth",2,"FontSize",12)
xlabel('t'), ylabel('x')
grid on
fprintf("Root-mean-squared error = %f\n", rmse(x,x_exac))