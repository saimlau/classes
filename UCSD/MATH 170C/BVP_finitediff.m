function [x,t]=BVP_finitediff(u,v,w,a,b,alpha,beta,h)
    t = (a:h:b)';
    n = length(t)-2;
    x = zeros(n,1); x0 = alpha; xe = beta;
    A = zeros(n,n);
    a = -1-h/2.*w(t); d = 2+h^2.*v(t); c = -1+h/2.*w(t);
    b = -1*h^2.*u(t(2:n+1)); b(1) = b(2)-a(1)*alpha; b(end) = b(end)-c(n+1)*beta;
    A(1,1:2) = [d(2) c(2)];
    for i=2:n-1
        A(i,(i-1):(i+1)) = [a(i) d(i+1) c(i+1)];
    end
    A(n,(n-1):n) = [a(n) d(n+1)];
    tem = A\b;
    x = [x0; tem; xe];
end