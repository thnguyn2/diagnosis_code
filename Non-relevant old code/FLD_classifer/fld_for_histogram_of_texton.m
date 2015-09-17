clc;
clear all;
close all;
load('H:\QPI data\texdir\33T19_texton_hist_40.mat','histim','lblim');
lblim1=lblim;
 load('H:\QPI data\label\33T19_resized.mat');
% lblim=lblim;
% figure(1);
% subplot(121);
% imagesc(lblim1);
% subplot(122);
% imagesc(lblim2);
% mask = ((lblim1.*lblim2)>0);
% %lblim =lblim1.*mask;
% lblim=lblim2;
ntextons = 51;
glandidx = find(lblim==1);
stromaidx = find(lblim==2);
lumenidx = find(lblim==0);
ngland = length(glandidx);
nstroma = length(stromaidx);
gland_hist = zeros(ngland,ntextons);
stroma_hist = zeros(nstroma,ntextons);
for bandidx=1:ntextons
    bandidx
    curband = histim(:,:,bandidx);
    glandvect=curband(glandidx);
    stromavect=curband(stromaidx);
    lumenvect=curband(lumenidx);
    gland_hist(:,bandidx)=glandvect(:);
    stroma_hist(:,bandidx)=stromavect(:);
end
clear histim;
clear curband;
clear glandvect;
clear stromavect;

gland_hist=gland_hist';
stroma_hist=stroma_hist';
% 
% %Load another image and combine
% load('H:\QPI data\texdir\44C19_texton_hist.mat','histim','lblim');
% ntextons = 80;
% glandidx2 = find(lblim==1);
% stromaidx2 = find(lblim==2);
% lumenidx2 = find(lblim==0);
% ngland2 = length(glandidx);
% nstroma2 = length(stromaidx);
% gland_hist2 = zeros(ngland,ntextons);
% stroma_hist2 = zeros(nstroma,ntextons);
% for bandidx=1:ntextons
%     bandidx
%     curband = histim(:,:,bandidx);
%     glandvect=curband(glandidx);
%     stromavect=curband(stromaidx);
%     lumenvect=curband(lumenidx);
%     gland_hist2(:,bandidx)=glandvect(:);
%     stroma_hist2(:,bandidx)=stromavect(:);
% end
% clear histim;
% 
% gland_hist2=gland_hist2';
% stroma_hist2=stroma_hist2';
% 
% gland_hist = [gland_hist gland_hist2];
% stroma_hist = [stroma_hist stroma_hist2];
% clear gland_hist2;
% clear stroma_hist2;
% 
% ngland = ngland + ngland2;
% nstroma = nstroma + nstroma2;

mugland = sum(gland_hist,2)/ngland;
mustroma = sum(stroma_hist,2)/nstroma;
mu = sum([gland_hist stroma_hist],2)/(ngland+nstroma);
figure(3)
plot(mugland)
hold on
plot(mustroma,'r')
Sglands = (gland_hist - repmat(mugland,[1 ngland]))*(gland_hist - repmat(mugland,[1 ngland]))';
Sstromas = (stroma_hist - repmat(mustroma,[1 nstroma]))*(stroma_hist - repmat(mustroma,[1 nstroma]))';
Sw=Sglands+Sstromas;
Sb = ngland*(mugland-mu)*(mugland-mu)'+nstroma*(mustroma-mu)*(mustroma-mu)';
[V,D]=eigs(Sb,Sw,2);

figure(5);
for sampleidx=1:100:300000
    plot(gland_hist(:,sampleidx),'b');
    hold on;
    plot(stroma_hist(:,sampleidx),'r');
end
legend('gland','stroma');
%Now, take only the first data and do the projections
nfeatures = 1;
W=real(V(:,1:nfeatures));
glandproj = W'*gland_hist;
stromaproj = W'*stroma_hist;

if (nfeatures==1)
    figure(4);
    hold on;
    plot(stromaproj','r');
     plot(glandproj','b');
   
end
legend('Stroma','Glands');
%Draw the projection map between grand and stroma


%Plot the distribution of two projections
if (nfeatures==2)
    figure(4);
    for glandidx=1:nglands
        plot(glandproj(1,glandidx),glandproj(2,glandidx),'+r');
        hold on;
    end
    for stromaidx=1:nstromas
        plot(stromaproj(1,stromaidx),stromaproj(2,stromaidx),'ob');
        hold on;
    end
    
end

%Draw the projection map between stroma & glands
clear gland_hist;
clear stroma_hist;
clear lumen_hist;

clear histim
load('H:\QPI data\texdir\33T19_texton_hist_40.mat','histim','lblim');
for textonidx=1:ntextons
    proj(:,:,textonidx) = histim(:,:,textonidx)*W(textonidx);
end

proj = sum(proj,3);
%Mask out the lumen
idx = find(lblim==0);
proj(idx)=0;
figure(1);
imagesc(abs(proj));
%imshow(abs(proj),[8e-3 18e-3])
colormap jet
colorbar;
title('1D reduction')
figure(2)
imagesc(lblim)
a=(proj);
glandidx = find(lblim==1);
stromaidx = find(lblim==2);
histbin = linspace(0,0.05,80);
glandhist = hist(a(glandidx),histbin)/length(glandidx);
stromahist = hist(a(stromaidx),histbin)/length(stromaidx);
figure(6);
plot(histbin,glandhist,'b');
hold on;
plot(histbin,stromahist,'r');
axis on;
legend('Gland','Stroma');
grid on;
hold off