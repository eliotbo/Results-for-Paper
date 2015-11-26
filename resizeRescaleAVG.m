function Z = resizeRescaleAVG(Z)

s = size(Z);
step1 = s(1)/256;
step2 = s(2)/256;

% Z = Z-min(min(Z));
M = imresize(Z,[256 256]);


Z = M;

Z = round(Z/max(max(Z))*255);