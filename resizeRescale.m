function Z = resizeRescale(Z)

s = size(Z);
step1 = s(1)/256;
step2 = s(2)/256;

% Z = Z-min(min(Z));
% Z = imresize(Z,[256 256]);

a=0;
for i=round(1:step1:s(1))
    a=a+1;aa=0;
    for ii=round(1:step2:s(2))
        aa=aa+1;
        M(a,aa) = Z(i,ii);
    end
end
Z = M;

Z = round(Z/max(max(Z))*255);