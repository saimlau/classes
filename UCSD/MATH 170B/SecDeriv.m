function ydd = SecDeriv(x,y)
h = x(2)-x(1);
n = length(x);
if n<4, error("More points are needed."), end
ydd = zeros(1,n);
ydd(1) = (-y(4)+4*y(3)-5*y(2)+2*y(1))/h^2;
for i=2:(n-1)
    ydd(i) = (y(i+1)-2*y(i)+y(i-1))/h^2;
end
ydd(n) = (-y(n-3)+4*y(n-2)-5*y(n-1)+2*y(n))/h^2;
end