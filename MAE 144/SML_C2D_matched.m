function Dz = SML_C2D_matched(Ds,h,omega,strictly_causal)
% function Dz = SML_C2D_matched(Ds,h,omega,strictly_causal)
% convert Ds(s) to Dz(z) using pole-zeromapping,
% h = timestep, omega = omega bar, strictly_causal = true/false
% Ds must be an RR_tf object
% TEST: SML_C2D_matched_test, it compares this code with matlab c2d(Ds,h,'matched'),
%       it also outputs symbolic results for comparsion with hand calculations
   if logical(sum(Ds.p==0)) && omega == 0, error("Omega cannot be 0 when there is a pole at 0!!!"), end
   if nargin<4, strictly_causal=false; if nargin<3, omega=0;end,if nargin<2,h=1;end,end
   syms temp
   m=Ds.num.n; n=Ds.den.n;
   c = string(class(Ds.z)); g = string(class(Ds.p));
   if c=="sym" || g=="sym", z_ = [temp]; p_ = [temp]; 
   else z_=[1]; p_=[1]; end
   if strictly_causal, z_(1:n-1) = ones(1,n-1)*-1; else, z_(1:n) = ones(1,n)*-1; end
   p_(1:n) = ones(1,n);
   for j = 1:n, p_(j) = exp(Ds.p(j)*h); end
   if m>=1, for j = 1:m, z_(j) = exp(Ds.z(j)*h);end, end
   tem = 1;
   for j=1:m, tem = tem*(omega*1i-Ds.z(j)); end
   for j=1:n, tem = tem/(omega*1i-Ds.p(j)); end
   tem2 = 1;
   for j=1:length(z_), tem2 = tem2*(exp(omega*1i*h)-z_(j)); end
   for j=1:length(p_), tem2 = tem2/(exp(omega*1i*h)-p_(j)); end
   K = norm(tem)/norm(tem2)*Ds.K;
   if isnan(K), error("Use a sllightly different omega bar!!!"), end
   if c=="sym" || g=="sym", z_=simplify(z_); p_=simplify(p_); K=simplify(K); end
   Dz = RR_tf(z_,p_,K);
   Dz.h = h;
end