clear,clc,close all
L = 3;
phi = 3*pi/5;
ns = 0:3:6;
xs = -2*pi:0.001:2*pi;
N = max(ns)+3;
for n=ns
    % f = @(x) cos(2*pi*n/L.*(x+phi));
    % plot(xs,f(xs))
    hold on
    % fh = @(x) 0;
    % Rh = @(x) 0;
    % for m=-N:N
    %     k = 2*pi/L*m;
    %     fh = @(x) fh(x)+fm(n,m,L,phi).*exp(1i*k.*x);
    %     Rh = @(x) Rh(x)+Rm(n,m).*exp(1i*k.*x);
    % end
    % plot(xs,Rh(xs),LineWidth=2)
    % plot(xs,fh(xs),"r--")
    RR = @(l) cos(2*pi*n/L*l);
    plot(xs,RR(xs),LineWidth=2)
end
xlim([-1,1])
ylim([-1.2,1.2])
ylabel("$R(\ell)$", Interpreter="latex")
xlabel("$\ell$", Interpreter="latex")
legend("n=0","n=3","n=6",Location="east")
set(gca,"LineWidth",2,"FontSize",15,"Box","on")

function out = fm(n,m,L,phi)
    out = 0;
    if n==m
        if n==0
            out = 1;
        else
            out = exp(1i*2*pi*n/L*phi)/2;
        end
    elseif n==-m
        out = exp(-1i*2*pi*n/L*phi)/2;
    end
end

function out = Rm(n,m)
    out = 0;
    if n==m
        if n==0
            out = 1;
        else
            out = 1/4;
        end
    end
end