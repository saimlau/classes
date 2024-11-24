clear,clc,close all
M = 27;
N = 512;
tsamp = ...
[    1;
    41;
    43;
    75;
    86;
   141;
   148;
   173;
   200;
   208;
   228;
   239;
   254;
   275;
   279;
   280;
   298;
   311;
   315;
   348;
   360;
   371;
   396;
   403;
   451;
   471;
   512];
ysamp = ...
[   0.0700;
    0.0274;
    0.0481;
   -0.1520;
   -0.2702;
   -0.0181;
   -0.1491;
   -0.0400;
   -0.1022;
   -0.1693;
   -0.0310;
    0.0631;
   -0.0194;
   -0.0307;
    0.0376;
    0.0532;
    0.3618;
    0.2472;
    0.1519;
   -0.0521;
    0.1546;
    0.1723;
   -0.0659;
   -0.0340;
   -0.0329;
    0.0806;
   -0.0252];

sigma = 0.003;

SHOWGRAPH = 1;
if SHOWGRAPH
	figure()
	plot(tsamp, ysamp, 'ko');
	axis([0 512 -0.5 0.5])
	legend('Sampled points')
end

%Problem 2
A = zeros([M,N])./sqrt(N);
for i = 2:N
    A(:,i) = sqrt(2/N).*cos(pi/2/N.*(2.*tsamp-1).*(i-1));
end
Yhat = A\ysamp;
for k=1:N
    tem = zeros([N,1]);
    tem(1:k) = Yhat(1:k);
    yhat = idct(tem,Type=2);
    if rmse(ysamp,yhat(tsamp))<3*sigma
        disp(k)
        break
    end
end

% To create a plot with yhat and ysamp together, uncomment the following code.
figure()
plot(tsamp, ysamp, 'ko');
hold on;
plot(1:N, yhat);
%% Mark the point yhat(129).
plot(129, yhat(129), 'rs')
