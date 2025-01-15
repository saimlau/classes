clear,clc,close all
NN = [21,51,101];
figure
set(gcf, 'Position',  [100, 100, 1500, 600])
xs = -pi:0.0001:pi;
f = @(x) x.*(x+2).*sin(x);
x0 = pi/5;
fx0 = f(x0)
for j=1:3
    N = NN(j);
    eval("P"+string(N)+"= @(x) sin(("+string(N)+"+1/2).*x)./sin(x./2);")
    eval("g"+string(N)+"= @(x) P"+string(N)+"(x-x0)./(2*pi).*f(x);")
    subplot(1,3,j)
    plot(xs,eval("P"+string(N)+"(xs)")./(2*pi))
    title("N="+string(N))
    xlim([-pi,pi])
    xticks(-pi:pi:pi)
    xticklabels(["-π","0","π"])
    set(gca,"FontSize",15,"LineWidth",2)
    eval("inte"+string(N)+"=abs(trapz(xs,g"+string(N)+"(xs))-fx0)")
end
