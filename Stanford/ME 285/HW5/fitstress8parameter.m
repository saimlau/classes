function sse = fitstress8parameter(c)
% <^>
%-----
%(~ ~)
% /(` 
%   )
%   Code to fit a biaxial experiment using following constituents:
%   --4 fiber collagen family (axial, circumferential, +/- 45deg),
%   --circumferential smooth muscle
%   -- elastin
%   ref paper -- Ramachandra, A.B., et al., Journal of biomechanical
%   engineering, 2015
%   questions - email - amarsden@stanford.edu
%   code for class purposes only
%   C index 
%   1-c1_mc,2-c2_mc
%   3-c1_collagen1,4-c2_collagen1
%   5-c1_collagen2,6-c2_collagen2
%   7-c1_collagen3,8-c2_collagen3
%   9-c1_collagen4,10-c2_collagen4
% 

 
lz=1.5748971193;
lt=1.0604937526:0.01:1.754257212;

%Fibre Angle
alpha_mp=pi/2;
alpha_co1=0;
alpha_co2=pi/2;
alpha_co3=c(8);
alpha_co4=-c(8);

%Prestretch
G_e_cir=1.4;
G_e_axi=1.6;
G_e=1.4;
G_c=1.08;
G_mp=1.2;

%mass fractions
 m_c=.45; m_e=0.12;m_m=.43;
 m_c1=.1*m_c;m_c2=.1*m_c;m_c3=.4*m_c;m_c4=.4*m_c;

% l - lambda, lz - lambda_z, lt - lambda_theta
l_mp=sqrt((lz*cos(alpha_mp)).^2+(lt*sin(alpha_mp)).^2);
l_co1=sqrt((lz*cos(alpha_co1)).^2+(lt*sin(alpha_co1)).^2);
l_co2=sqrt((lz*cos(alpha_co2)).^2+(lt*sin(alpha_co2)).^2);
l_co3=sqrt((lz*cos(alpha_co3)).^2+(lt*sin(alpha_co3)).^2);
l_co4=sqrt((lz*cos(alpha_co4)).^2+(lt*sin(alpha_co4)).^2);
%  -----------------------------------------Thetha Stress------------------
%Passive SMC
W_mp=c(1)*(exp(c(2)*(G_mp^2*l_mp.^2-1).^2));
dexpdlt_mp=c(2)*(G_mp^4*(4*(lt.^3).*(sin(alpha_mp)).^4+4*lt.*(lz.^2).*(sin(alpha_mp)).^2*(cos(alpha_mp)).^2)-G_mp^2*(4*lt.*(sin(alpha_mp)).^2));

dWdlt_mp=W_mp.*dexpdlt_mp;%

%Collagen 1 at 0 degrees
W_co1=c(3)*(exp(c(4)*(G_c^2*l_co1.^2-1).^2));
dexpdlt_co1=c(4)*(G_c^4*(4*(lt.^3).*(sin(alpha_co1)).^4+4*lt.*lz.^2*(sin(alpha_co1)).^2*(cos(alpha_co1)).^2)-G_c^2*(4*lt.*(sin(alpha_co1)).^2));

dWdlt_co1=W_co1.*dexpdlt_co1;

%Collagen 2 
W_co2=c(1)*(exp(c(2)*(G_c^2*l_co2.^2-1).^2));
dexpdlt_co2=c(2)*(G_c^4*(4*lt.^3*(sin(alpha_co2)).^4+4*lt.*lz.^2*(sin(alpha_co2)).^2*(cos(alpha_co2)).^2)-G_c^2*(4*lt.*(sin(alpha_co2)).^2));

dWdlt_co2=W_co2.*dexpdlt_co2;

%Collagen 3
W_co3=c(5)*(exp(c(6)*(G_c^2*l_co3.^2-1).^2));
dexpdlt_co3=c(6)*(G_c^4*(4*(lt.^3).*(sin(alpha_co3)).^4+4*lt.*lz.^2*(sin(alpha_co3)).^2*(cos(alpha_co3)).^2)-G_c^2*(4*lt.*(sin(alpha_co3)).^2));

dWdlt_co3=W_co3.*dexpdlt_co3;

%Collagen 4
W_co4=c(5)*(exp(c(6)*(G_c^2*l_co4.^2-1).^2));
dexpdlt_co4=c(6)*(G_c^4*(4*lt.^3*(sin(alpha_co4)).^4+4*lt.*lz.^2*(sin(alpha_co4)).^2*(cos(alpha_co4)).^2)-G_c.^2*(4*lt*(sin(alpha_co4)).^2));

dWdlt_co4=W_co4.*dexpdlt_co4;

% Elastin
dWdlt_e=c(7)*(2.*lt.*G_e_cir^2-2./(G_e_cir^4*lt.^3.*lz.^2));


%s_theta=lt.*(m_e*dWdlt_e+m_c1*dWdlt_co1+m_c2*dWdlt_co2+m_c3*dWdlt_co3+m_c4*dWdlt_co4+m_m*dWdlt_mp);
s_theta=lt.*(dWdlt_e+dWdlt_co1+dWdlt_co2+dWdlt_co3+dWdlt_co4+dWdlt_mp);


%------------------------------------ Z Stress-----------------------------

%Passive SMC
W_mp=c(1)*(exp(c(2)*(G_mp^2*l_mp.^2-1).^2));
dexpdlz_mp=c(2)*(G_mp^4*(4*(lz.^3).*(cos(alpha_mp)).^4+4*lz.*(lt.^2).*(sin(alpha_mp)).^2*(cos(alpha_mp)).^2)-G_mp^2*(4*lz.*(cos(alpha_mp)).^2));

dWdlz_mp=W_mp.*dexpdlz_mp;%

%Collagen 1 at 0 degrees
W_co1=c(3)*(exp(c(4)*(G_c^2*l_co1.^2-1).^2));
dexpdlz_co1=c(4)*(G_c^4*(4*(lz.^3).*(cos(alpha_co1)).^4+4*lz.*lt.^2*(sin(alpha_co1)).^2*(cos(alpha_co1)).^2)-G_c^2*(4*lz.*(cos(alpha_co1)).^2));

dWdlz_co1=W_co1.*dexpdlz_co1;

%Collagen 2 
W_co2=c(1)*(exp(c(2)*(G_c^2*l_co2.^2-1).^2));
dexpdlz_co2=c(2)*(G_c^4*(4*(lz.^3)*(cos(alpha_co2)).^4+4*lz.*lt.^2*(sin(alpha_co2)).^2*(cos(alpha_co2)).^2)-G_c^2*(4*lz.*(cos(alpha_co2)).^2));

dWdlz_co2=W_co2.*dexpdlz_co2;

%Collagen 3
W_co3=c(5)*(exp(c(6)*(G_c^2*l_co3.^2-1).^2));
dexpdlz_co3=c(6)*(G_c^4*(4*(lz.^3).*(cos(alpha_co3)).^4+4*lz.*lt.^2*(sin(alpha_co3)).^2*(cos(alpha_co3)).^2)-G_c^2*(4*lz.*(cos(alpha_co3)).^2));

dWdlz_co3=W_co3.*dexpdlz_co3;

%Collagen 4
W_co4=c(5)*(exp(c(6)*(G_c^2*l_co4.^2-1).^2));
dexpdlz_co4=c(6)*(G_c^4*(4*(lz.^3).*(cos(alpha_co4)).^4+4*lz.*lt.^2*(sin(alpha_co4)).^2*(cos(alpha_co4)).^2)-G_c^2*(4*lz.*(cos(alpha_co4)).^2));

dWdlz_co4=W_co4.*dexpdlz_co4;

% Elastin

dWdlz_e=c(7)*(2.*lz.*G_e_axi^2-2./(G_e_axi^4*lt.^2.*lz.^3));

%s_z=1./lt.*(m_e*dWdlz_e+m_c1*dWdlz_co1+m_c2*dWdlz_co2+m_c3*dWdlz_co3+m_c4*dWdlz_co4+m_m*dWdlz_mp);
s_z=lz.*(dWdlz_e+dWdlz_co1+dWdlz_co2+dWdlz_co3+dWdlz_co4+dWdlz_mp);

%Plot results
figure(1);plot(lt,s_theta,'LineWidth',2);
xlabel('\lambda_\theta');ylabel('\sigma_\theta(Pa)')
hold on
figure(2);plot(lt,s_z,'LineWidth',2);
xlabel('\lambda_\theta');ylabel('\sigma_z(Pa)')
hold on

end
% <^>
%-----
%(~ ~)
% /(` 
%   )
