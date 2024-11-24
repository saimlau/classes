function [x,t] = AM4(f,x0,a,b,h,TOL,MaxIters)
    t = a:h:b;
    n = length(t);
    x = zeros(1,n);
    [x(1:3), ~] = RK4(f,x0,t(1),t(3),h);
    for i = 3:n-1
        xi = x(i); ti = t(i);
        g = @(y) xi+h/24*(9*f(t(i+1),y)+19*f(ti,xi)-5*f(ti-h,x(i-1))+f(ti-2*h,x(i-2)));
        y = xi;
        for j = 1:MaxIters
            if abs(y-g(y))>TOL
                y = g(y);
            else, break
            end
        end
        x(i+1) = y;
    end
end