%Optimization in the pixel basis (PB) of an image affected by Poisson
%noise.

%inputs: im, opt and plotsOn.
%im.name: name of the image in the folder
%im.scale: modifies size of image
%im.photons: number of photons per pixel at maximum brightness
%alpha:
%gam: control over convergence rate (0.1 is typical)
%lambda: control over sparsity of DCB representation
%beta: control over total variation
%huberSmoothing: set to 1 to use smoothing on absolute and 0 not to
%mu: Huber smoothing constant for sign of image in DCB representation
%mu2: Huber smoothing constant for sign of variations in pixel basis
%squareTotalVar: set to 1 to optimize square instead of absolute value
%background: number of photons in background (3 is more stable than 0)
%plotsOn: 1 will plot initial, noisy and recovered images and cost function
function [results M]= f1_pbBasis(N,Z,beta,squareTotalVar)

plotsOn = 0;

% im.name = 'cameraman.tif';
% im.scale = 1;
% im.photons = 20;

% beta = 0.5;
% squareTotalVar = 1;

alpha = 0.8;      
gam = 0.1;
lambda = 0.;
huberSmoothing=0;
mu = 1;
mu2 = 1;
numIt= 20;
background = 3;

s=size(Z);
r = s(1) ;
c = s(2) ;

%Add background to noisy image
 N = N+background;

%Discrete Cosine Transform matrices
Ur = dctmatrix(r);
Uc = dctmatrix(c);

% initial guess
M = N-background;

v = 0*M;
bestf=1E20;
lfacN =  log(factorial(N));
tic
for i=1:numIt
    %Anticipate that image will approximately go towards momentum direction
    M = M + alpha*v;
    
    %gradient of likelihood function in pixel basis
    gradientLikelihood =  (1-N./(M+background));
     fLikelihood(i) = sum(sum(-N.*log(M+background)+M+background + lfacN));
     
    %Calculate gradient of DCB sparsity cost function if lambda>0
    if lambda>0
        %basis change
        A = Ur*M*Uc';
        squaredTwoNorm = sum(sum(A.^2));
        oneNorm = sum(sum(abs(A)));
        
        %choose smoothing or not for abs(A)
        if huberSmoothing == 1
            signAA = huberAbs(A,mu);
        else
            signAA = sign(A);
        end
        
        oneTwo = oneNorm/squaredTwoNorm;
        gradientDCT = 2*oneTwo*(Ur'*signAA*Uc-oneTwo*M);
        fDCT(i) = lambda*oneNorm^2/squaredTwoNorm;
    else
        gradientDCT = 0;
        fDCT(i) = 0;
    end
    
    
    %Calculate gradient of total variation cost function if beta>0
    if beta>0
        Imx1 = circshift(M,[1 0]);
        Imy1 = circshift(M,[0 1]);
        Imx2 = circshift(M,[-1 0]);
        Imy2 = circshift(M,[0 -1]);
        
        subX = Imx1-M;
        subY = Imy1-M;
        subX2 = M-Imx2;
        subY2 = M-Imy2;
        
        %difference between first column (row) and last column (row) is
        %meaningless
        subX(end,:) = 0;
        subY(:,end) = 0;
        subX2(1,:) =  0;
        subY2(:,1) =  0;
        
        if squareTotalVar == 1
            ImSubx = -(subX);
            ImSuby = -(subY);
            ImSubx2 = (subX2);
            ImSuby2 = (subY2);
        else
            if huberSmoothing == 1
                ImSubx = -huberAbs(subX,mu2);
                ImSuby = -huberAbs(subY,mu2);
                ImSubx2 = huberAbs(subX2,mu2);
                ImSuby2 = huberAbs(subY2,mu2);
            else
                ImSubx = -sign(subX);
                ImSuby = -sign(subY);
                ImSubx2 = sign(subX2);
                ImSuby2 = sign(subY2);
            end
        end
        gradientTotalVar = ImSubx + ImSuby + ImSubx2 + ImSuby2;
        fTotalVar(i)  = beta*(sum(sum(abs(subX)))+sum(sum(abs(subY))));
    else
        gradientTotalVar = 0;
        fTotalVar(i)  = 0;
    end
    
    %total gradient
    gradient = gradientLikelihood+lambda*gradientDCT+beta*gradientTotalVar;
    
    %evaluate cost function
    f(i) = fLikelihood(i) + fDCT(i) + fTotalVar(i);
    
    %accelerated gradient descent
    M  = M - alpha*v;
    v = alpha*v - gam*gradient;
    M = M + v;
    if f(i) < bestf
        bestf = f(i);
        bestM = M;
    end
end

%If the problem is not convex, recover image with lowest cost function
M = bestM;
% M(M<0) = 0;
% M(M>2.5) = 2.5;

runningTime = toc;

PSNRf.noisy = PSNR(N-background,Z);
PSNRf.recovered = PSNR(M,Z);
totsPhotons = sum(sum(N-background));

results.totsPhotons = totsPhotons;
results.time = runningTime;
results.in = PSNRf.noisy;
results.out = PSNRf.recovered;


%%%%%%%%%%%% Plot final image and cost functions %%%%%%%%%%%%
if plotsOn == 1
figure(21)
imagesc(M); colormap gray

figure(11)
subplot(4,1,1); hold on
plot(real(f));  ylabel('f')

subplot(4,1,2);     hold on
plot(real(fDCT));   ylabel('fDCT')

subplot(4,1,3);             hold on
plot(real(fLikelihood));    ylabel('fLikelihood')

subplot(4,1,4);         hold on
plot(real(fTotalVar));  ylabel('Total Variance')
end
%%%%%%%%%%%% Plot final image and cost functions %%%%%%%%%%%%

