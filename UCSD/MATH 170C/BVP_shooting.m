function [x,t]=BVP_shooting(f,a,b,alpha,beta,h,TOL,MaxIters)
    g = @(tt,y) [y(2); f(tt,y(1),y(2))];
    z0 = 0; z1 = 1;
    le = phi(z0,g,a,b,alpha,beta,h);
    err = phi(z1,g,a,b,alpha,beta,h);
    for i=1:MaxIters
        if abs(err) >= TOL
            z2 = z1-err*(z1-z0)/(err-le);
            z0 = z1;
            z1 = z2;
            le = err;
            err = phi(z1,g,a,b,alpha,beta,h);
        else
            break
        end
    end
    [tem3,t] = RK4(g, [alpha; z1], a, b, h);
    x = tem3(1,:);
end

function k = phi(z,g,a,b,alpha,beta,h)
    [tem,~] = RK4(g, [alpha; z], a, b, h);
    k = tem(1,end)-beta;
end