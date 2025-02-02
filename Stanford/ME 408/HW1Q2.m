clear,clc,close all
N = 4;
kk = -N:N;
fh = @(x) 0;
gh = @(x) 0;
for k=-N:N
    fh = @(x) fh(x)+fk(k).*exp(1i*k.*x);
    gh = @(x) gh(x)+gk(k).*exp(1i*k.*x);
end
xs = -2*pi:0.01:2*pi;
plot(xs, fh(xs))
xlabel("x")
xticks(-2*pi:pi:2*pi)
xticklabels(["-2π","-π","0","π","2π"])
set(gca,"FontSize",15,"LineWidth",2)
figure
plot(xs, gh(xs))
xlabel("x")
xticks(-2*pi:pi:2*pi)
xticklabels(["-2π","-π","0","π","2π"])
set(gca,"FontSize",15,"LineWidth",2)

function out = fk(k)
% if mod(k,2)==0
%     out = -k/(k^3-k)/pi;
% else
%     out = 0;
% end
if abs(k)==1
    out = 1i*sign(k)/(-2*2);
else
    out = (cos(pi*k)+1)/(2*pi*(1-k^2));
end
end

function out = gk(k)
if k ~= 0
    out = 1i/k*(-1)^k;
else
    out = 0;
end
end