%This program demonstrate the dimensionality reduction for the textondata
clc;
clear all;
close all;
load texton_data.mat;
hist_data=hist_data_over_textons(:,1:3);
threshval = 0.6;
%Choose the threshold 
idx = find(max(hist_data,[],2)>threshval);
hist_data= hist_data(idx,:);
texton_map_over_images=texton_map_over_images(idx,:);
texton_map_over_images=repmat(sqrt(sum(texton_map_over_images.^2,2)),1,20).*texton_map_over_images;
figure
imagesc(hist_data);
[maxval,labelidx]=max(hist_data,[],2);
lumenidx = find(labelidx==1);
glandidx = find(labelidx==2);
stromaidx = find(labelidx==3);
figure(1);
for spectrumidx=1:length(glandidx)
    plot(texton_map_over_images(glandidx(spectrumidx),:));
    hold on;
end
title('Gland');
figure(2);
for spectrumidx=1:length(stromaidx)
    plot(texton_map_over_images(stromaidx(spectrumidx),:),'r');
    hold on;
end
title('Stroma');
% figure(3);
% for spectrumidx=1:length(lumenidx)
%     plot(texton_map_over_images(lumenidx(spectrumidx),:),'r');
%     hold on;
% end
title('Lumen')
%Perform FLD on the dataset to pickout the difference....
glanddata = texton_map_over_images(glandidx,:);
lumendata = texton_map_over_images(lumenidx,:);
stromadata = texton_map_over_images(stromaidx,:);
glanddata = glanddata';
lumendata =lumendata';
stromadata = stromadata';
mugl = mean(glanddata,2);
mustr = mean(stromadata,2);
mulum = mean(lumendata,2);
ngl = size(glanddata,2);
nstr = size(stromadata,2);
nlum = size(lumendata,2);
mu = (ngl*mugl + nstr*mustr + nlum * mulum)/(ngl+nstr+nlum);
sgl = glanddata*glanddata';
sstr = stromadata*stromadata';
slum = lumendata*lumendata';
%Compute the covariance matrix for within class variation
sw = sgl + sstr + slum;
%Between class scatterer
sb = ngl*(mugl-mu)*(mugl-mu)'+nstr*(mustr-mu)*(mustr-mu)'+...
    nlum*(mulum-mu)*(mulum-mu)';
[W,D]=eig(sb,sw);
%Now, do the projection on the first 4 vectors
nvect=1;
Wopt = W(:,1:nvect);
newgldata=Wopt'*glanddata;
newstrdata=Wopt'*stromadata;
newlumdata=Wopt'*lumendata;
figure(4);
for spectrumidx=1:length(glandidx)
    plot(newgldata(:,spectrumidx));
    hold on;
end
title('Gland');
figure(5);
for spectrumidx=1:length(stromaidx)
    plot(newstrdata(:,spectrumidx));
    hold on;
end
title('Stroma');



