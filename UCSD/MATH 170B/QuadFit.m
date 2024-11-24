function [a,Q] = QuadFit(x,y)
    syms xx
    k=2;
    n = length(x)-1;
    z = zeros(1,n);
    D = zeros(1,n-1);
    a = zeros(n-1,k+1);
    
    for i=1:n
        D(i) = y(i)-(x(i)-x(i+1))/2*z(i);
        z(i+1) = 2*(y(i+1)-D(i))/(x(i+1)-x(i));
        tem = z(i)/(2*(x(i)-x(i+1)))*(xx-x(i+1))^2+z(i+1)/(2*(x(i+1)-x(i)))*(xx-x(i))^2+D(i);
        if i==1
            Q = piecewise(xx<=x(i+1),tem);
        elseif i<n
            Q = piecewise(xx<x(i),Q,xx<=x(i+1),tem);
        else
            Q = piecewise(xx<x(i),Q,x(i)<=xx,tem);
        end
        atem = sym2poly(tem);
        a(i,(k+1-length(atem)+1):(k+1)) = sym2poly(tem);
    end
end