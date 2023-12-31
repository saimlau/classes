% 2a. Consider the (pathalogical...) plant G(s)=[(s+2)(s-2)(s+5)(s-5)]/[(s+1)(s-1)(s+3)(s-3)(s+6)(s-6)]=b(s)/a(s).  
% Compute a controller D(s)=y(s)/x(s) that puts the 6 poles of T(s)=G(s)D(s)/[1+G(s)D(s)]=g(s)/f(s) at 
% s={-1,-1,-3,-3,-6,-6}.  Check your answer to make sure it works.  Hint: form the denominator polynomial f(s), 
% then call RR_Diophantine.  The code that you need to write for this problem is literally a 2-liner: one to call 
% RR_Diophantine properly, then one to double check that it worked.  Note that, amongst all solutions of 
% a(s)*x(s)+b(s)*y(s)=f(s), RR_Diophantine returns the answer with the smallest order for y(s).
clear,clc,close all
a = RR_poly([-1 1 -3 3 -6 6],1); b = RR_poly([-2 2 -5 5],1); f = RR_poly([-1,-1,-3,-3,-6,-6],1); G = RR_tf(b,a);
[x1,y1] = RR_diophantine(a,b,f); D = RR_tf(y1,x1)
test = trim(a*x1+b*y1), residual = norm(test-f);