function plotSub(M,num,clims,r)

subplot(1,5,num);
imagesc(M,clims)
colormap gray
axis image
set(gca,'YTickLabel',[],'XTickLabel',[])
set(gca,'YTick',[],'XTick',[])

if num==1
    
elseif num==2
    title(['PNRS=' num2str(r.in,4) ' dB, ' num2str(r.totsPhotons,5) ' photons']);
else
    
title(['PNRS=' num2str(r.out,4) ' dB, ' num2str(r.time,3) 's']);
end