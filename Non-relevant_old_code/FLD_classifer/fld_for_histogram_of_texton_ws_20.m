clc;
clear all;
close all;
load('H:\QPI data\texdir\33T19_texton_hist_ws_20.mat','histim','lblim');
ntextons = 20;
glandidx = find(lblim==1);
stromaidx = find(lblim==2);
lumenidx = find(lblim==0);
ngland = length(glandidx);
nstroma = length(stromaidx);
nlumen = length(lumenidx);
gland_hist = zeros(ngland,ntextons);
stroma_hist = zeros(nstroma,ntextons);
lumen_hist = zeros(nlumen,ntextons);
for bandidx=1:ntextons
    bandidx
    curband = histim(:,:,bandidx);
    glandvect=curband(glandidx);
    stromavect=curband(stromaidx);
    lumenvect=curband(lumenidx);
    gland_hist(:,bandidx)=glandvect(:);
    stroma_hist(:,bandidx)=stromavect(:);
    lumen_hist(:,bandidx)=lumenvect(:);
end
clear histim;

gland_hist=gland_hist';
lumen_hist=lumen_hist';
stroma_hist=stroma_hist';
mugland = sum(gland_hist,2)/ngland;
mustroma = sum(stroma_hist,2)/nstroma;
mulumen = sum(lumen_hist,2)/nlumen;
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

%Now, take only the first data and do the projections
nfeatures = 1;
W=real(V(:,1:nfeatures));
glandproj = W'*gland_hist;
stromaproj = W'*stroma_hist;

if (nfeatures==1)
    figure(4);
    hold on;
    plot(stromaproj','r');
     plot(glandproj');
   
end
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
load('H:\QPI data\texdir\33AA12_texton_hist_ws_20.mat','histim','lblim');
for textonidx=1:ntextons
    histim(:,:,textonidx) = histim(:,:,textonidx)*W(textonidx);
end

proj = sum(histim,3);
%Mask out the lumen
idx = find(lblim==0);
proj(idx)=0;
figure(1);
imagesc(abs(proj));
colorbar;
title('1D reduction')
figure(2)
imagesc(lblim)

a=abs(proj);
glandidx = find(lblim==1);
stromaidx = find(lblim==2);
histbin = linspace(0,0.013,80);
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