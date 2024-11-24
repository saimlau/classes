% The following data were measured in a pig's main pulmonary artery.
clear,clc,close all
flow=[33.42
    56.19
    73.697
    96.721
    139.85
    164.46
    177.44
    196.25
    198.77
    184.72
    162.09
    131.85
    91.057
    75.404
    62.991
    32.539
    21.865
    28.182
    23.896
    19.457
    19.911
    13.432
    5.284
    -1.0584]; % cm^3 / second

period=0.594; % seconds

radius=1.32; % centimeters
rho = 1.06; %g/cm^3
nu = 0.04; %dyne*s/cm^2

N = length(flow);
omega = 2*pi/period;
Q = fft(flow);
ts = 0:period/8:period;
rs = -1:0.01:1;
rrs = rs.*radius;
nr = length(rs);
nt = length(ts);

alpha = @(n) radius.*sqrt(n.*omega./nu);
tem1 = @(r,n) (besselj(0,alpha(n).*1i^(3/2))-besselj(0,alpha(n).*r./radius.*1i^(3/2)))./(alpha(n).*1i^(3/2).*besselj(0,alpha(n).*1i^(3/2))-radius.*besselj(1,alpha(n).*1i^(3/2)));
tem2 = @(r,n,t) B(n,N,Q).*alpha(n).*1i.^(3/2).*tem1(r,n).*exp(1i.*n.*omega.*t);
v = @(r,t) real(sum(tem2(r,[-N/2:-1 1:N/2],t)))/(2*pi*radius);
vs = zeros(nt,nr);
for it=1:nt
    for ir=1:nr
        vs(it,ir) = v(rrs(ir),ts(it));
    end
end
figure()
for it=1:nt
    subplot(3,3,it)
    plot(vs(it,:),rs)
    title("t="+ts(it))
    xlabel("velocity (cm/s)")
    ylabel("r/R")
    xlim([-10 15])
end
fprintf("w(0,0) = %.3f cm/s\n", v(0,0))

function Bn = B(n,N,Q)
    Bn = zeros(1,length(n));
    for j=1:length(n)
        if abs(n(j))==N/2
            Bn(j) = Q(N/2+1)/2/N;
        elseif n(j)<0
            Bn(j) = Q(n(j)+N+1)/N;
        else
            Bn(j) = Q(n(j)+1)/N;
        end
    end
end
