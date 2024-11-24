function [x,t] = RK4(f,x0,a,b,h)
    t = a:h:b;
    n = length(t);
    m = length(x0);
    x = zeros(m,n); x(:,1) = x0;
    for i=1:n-1
        ti = t(i);
        xi = x(:,i);
        X1 = xi;
        X2 = xi+h/2.*f(ti,X1);
        X3 = xi+h/2.*f(ti+h/2,X2);
        X4 = xi+h.*f(ti+h/2,X3);
        x(:,i+1) = xi+h/6.*(f(ti,X1)+2.*f(ti+h/2,X2)+2.*f(ti+h/2,X3)+f(ti+h,X4));
    end
end