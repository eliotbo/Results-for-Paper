clear;clc;

addpath('./PURE-LET');


name{1}  = 'cameraman.tif';
name{2} = 'barbara.png';
name{3} = 'moon.gif';

BETA =   [0.4 0.4 0.7 ];
LAMBDA = [0.003 0.0015 0.002 ];
SQUARE = [0 0 1 ];
SCALE = [0.5 1 1];
PHOTONS = [10 10 5 ];

images = [];
for pos = [1:3]

    beta = BETA(pos);
    lambda = LAMBDA(pos);
    square = SQUARE(pos);
    scale = SCALE(pos);
    
nameIm = name{pos};

photonPerPixel = PHOTONS(pos);

Z = im2double(imresize(imread(nameIm),scale));
if pos==3
    Z = Z(1:700,1:700);
end
Z = Z/max(max(Z));
Z = Z*photonPerPixel;

N = poissrnd(Z);

[PIXresults PIXim]= f1_pbBasis(N,Z,beta,square);
[DCTresults DCTim]= f2_dctShrink(N,Z,lambda);
[LETresults LETim]= purelet_denoising_1(Z,N);

LETresults.out = PSNR(LETim,Z);

M{1,pos} = Z;
M{2,pos} = N;
M{3,pos} = LETim;
M{4,pos} = PIXim;
M{5,pos} = DCTim;

metric(pos,1:4) = [PIXresults.in LETresults.out PIXresults.out DCTresults.out];
times(pos,1:3) = [LETresults.time PIXresults.time DCTresults.time];
numPhots(pos) = PIXresults.totsPhotons;

bufferV = ones(256,3)*255;


MM=[];
for i=1:5
    if i==1
         M{i,pos} = resizeRescaleAVG(M{i,pos});
    else
        
        M{i,pos} = resizeRescale(M{i,pos});
    end
    MM = [MM  bufferV M{i,pos} bufferV];
    
    
end

s= size(MM);
sizeh = s(2);
bufferH = ones(3,sizeh)*255;

images = [images; bufferH; MM ; bufferH];

end

% clims = [0 photonPerPixel+sqrt(photonPerPixel)/4];
% 
% figure(21)
% plotSub(Z,1,clims,PIXresults)
% plotSub(N,2,clims,PIXresults)
% plotSub(LETim,3,clims,LETresults)
% plotSub(PIXim,4,clims,PIXresults)
% plotSub(DCTim,5,clims,DCTresults)

figure(65)
imagesc(images)

colormap gray
axis image
set(gca,'YTickLabel',[],'XTickLabel',[])
set(gca,'YTick',[],'XTick',[])

metric
times
numPhots
% matrix2latex([1 2; 3 4],'ge')
