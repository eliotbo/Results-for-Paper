function [results M]= f2_dctBasis(N,A)

Z =A;
s=size(A);
d=s(1);
d2=s(2);

background = 0.;
Anoisy = N+background;

D = dctmtx(size(A,1));
D2 =dctmtx(size(A,2)); 

d =length(Anoisy);

%make matrix into a vector
u = reshape(Anoisy,[prod(size(A)) 1]);

%plus and minus basis change matrices
Dp = D.*(D>=0);
Dm = D.*(D<0);
Tp = Dp*D';
Tm = Dm*D';

tic
%DCT of data
N =D*Anoisy*D2';

% initial guess
M = N;

S = abs(D)*Anoisy*abs(D2)';



alpha = 0.7;
gam = 10;
lambda1 = 0.07;
calculateCost = 1;
huber = 0;
 mu = 0.5;

a = 0;
numIt=  50;
v = 0*M;

for i=1:numIt
    gradientLikelihood =   N./(S+M) - M./(N+1/2+sqrt(S.^2-M.^2+(N+3/2).^2));

    twoNorm = sum(sum(abs(M).^2));
    squaredOneNorm = sum(sum(abs(M)))^2;
    
    signM = sign(M);
    if huber == 1
    signM  = signM.*(abs(M)>mu);
    signM =signM +  (abs(M)<mu).*M/mu;
    end
    
    grad1Norm = 2*sum(sum(abs(M))).*signM;
    grad2Norm = M;

    gradientDCT = grad1Norm/twoNorm - squaredOneNorm/twoNorm^2*grad2Norm;

%     gradientDCT = signM;
    gradient = -gradientLikelihood + lambda1*gradientDCT;

    fDCT(i) = lambda1*squaredOneNorm/twoNorm;
    if calculateCost == 1 && mod(i,2)==0
        a=a+1;
        fLikelihood(a) = sum(sum( N/2.*log((S+M+eps)./(S-M+eps))+log(1E-9+besseli(N,sqrt(S.^2-M.^2)))));
        f(a) = fLikelihood(a) + fDCT(a);
    end

    v = alpha*v - gam*gradient;
    M = M + v;

end
results.time = toc;


finalI = D'*M*D2;
finalI(finalI<0) = 0;

M = finalI;

results.in = PSNR(Anoisy,A);
results.out = PSNR(finalI,A);





