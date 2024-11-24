clear,clc,close all
x = [-7,-4,-1,0,2,5,7];
y = [20,14,5,3,-2,-10,-15];
[a,~] = QuadFit(x,y)
plotQuadFit(x,y,a)


function plotQuadFit(x,y,a)
    syms xx
    n = length(x)-1;
    for i=1:n
        tem = poly2sym(a(i,:),xx);
        if i==1
            Q = piecewise(xx<=x(i+1),tem);
        elseif i<n
            Q = piecewise(xx<x(i),Q,xx<=x(i+1),tem);
        else
            Q = piecewise(xx<x(i),Q,x(i)<=xx,tem);
        end
    end
    figure()
    fplot(Q,[x(1)-2, x(end)+2],'g-', 'LineWidth',3), hold on
    plot(x,y,'r.','MarkerSize',15)
    legend('Quad spline $$Q(x)$$','$${(x_i,y_i)}_{i=0}^n$$','Interpreter','latex')
    set(gca,"FontSize",15,"LineWidth",1.5,'TickLabelInterpreter','latex')
    xlabel('x','Interpreter','latex')
    ylabel('y','Interpreter','latex')
    axis equal
    ylim([min(y)-1,max(y)+1])
end

