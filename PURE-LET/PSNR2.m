function result = PSNR(M1,M2)

s = size(M1);
% MAXi = max(max([M1 M2]));
MAXi = 255;
% M1 = M1*MAXi/max(max(M1));
% M2 = M2*MAXi/max(max(M2));

m1 = mean(mean(M1));
m2 = mean(mean(M2));

M1 = m2/m1*M1;

MSE = 1/prod(s)*sum(sum((M1-M2).^2));

% MSE = 1/prod(s)*norm((M1-M2),'fro')^2;

result = 10*log(MAXi^2/MSE);