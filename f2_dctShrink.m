function [results M]= f2_dctShrink(N,Z,lambda)

s=size(N);
r=s(1);
c=s(2);
photonPerPixel = max(max(N));

Ur = dctmtx(r);
Uc =dctmtx(c); 

% lambda = 0.003;

tic;
S = abs(Ur)*N*abs(Uc)';
T =Ur*N*Uc';
A1 = T - S.*lambda.*sign(T);
sameSign = (1+sign(A1).*sign(T))/2;
%if element of Te has a different sign from T, set element of Te to 0.
A2 = A1.*sameSign;
M = Ur'*A2*Uc;
results.time = toc;

results.in = PSNR(N,Z);
results.out = PSNR(M,Z);







