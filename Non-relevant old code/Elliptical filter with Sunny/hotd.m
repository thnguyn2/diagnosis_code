%This program demonstrates the use of ellipse gaussian kernel to generate
%the histogram of gradient response from the edir images
clear all;
clc;
close all;
load edir;
im = edir(:,:,5);
ndir = 8;
Fsym=makefilterfortexdir(1.0);
Fnonsym= makefilterfortexdir(3.0);
%Generate
figure
imagesc(im);
nrows =size(im,1);
ncols =size(im,2);
hist = zeros(nrows,ncols,ndir);
histfiltervalsym = zeros(nrows,ncols,ndir);


%Form a new image that has coherent values of the texture direction accross
%scale
nscales = 5;
hist = zeros(nrows,ncols,ndir);
for idx = 1:ndir
    finhist = ones(nrows,ncols);
    for scaleidx=1:nscales
        mask = zeros(nrows,ncols);
        curim = edir(:,:,scaleidx);
        pixidx = find(curim==idx); %Find satisfied pixels at current scale
        mask(pixidx)=1;
        finhist = finhist.*mask;
    end
    hist(:,:,idx)=finhist;
end
goodmask=sum(hist,3);
newedir = edir(:,:,5).*goodmask;
figure(1);
imagesc(newedir);
drawnow;
for idx = 1:ndir
    histfiltervalsym(:,:,idx)=imfilter(hist(:,:,idx),Fsym(:,:,idx),'same'); 
end

%Find the index of the best direction for the window
[maxval,maxidx]=max(histfiltervalsym,[],3);
allpixelidx=1:nrows*ncols;

%Do an extra layer of boosting by filtering the surrounding
%Now, compute the histogram over the new window
histfiltervalnonsym = zeros(nrows,ncols,ndir);
for idx = 1:ndir
    idx
     pixidx = find(newedir ==idx);
     mask = zeros(nrows,ncols);
     mask(pixidx)=1;
     %Filer each image with all posible direction
     tempcube = zeros(nrows,ncols,ndir);
     for diridx = 1:ndir
        tempcube(:,:,diridx) = imfilter(mask,Fnonsym(:,:,diridx),'same');
     end
     %Pick out the response corresponding to the best window
     bestidx = (maxidx(allpixelidx)-1)*nrows*ncols+allpixelidx;%Compute the 3D index
     bestim = zeros(nrows,ncols);
     bestim(allpixelidx)=tempcube(bestidx);
     histfiltervalnonsym(:,:,idx)=bestim;
end
%Normalization to a probability vector
sumim =1./sum(histfiltervalnonsym,3);
histfiltervalnonsym = histfiltervalnonsym.*(repmat(sumim,[1 1 ndir]));
[maxvalnonsym,maxidxnonsym]=max(histfiltervalnonsym,[],3);
figure(2);
subplot(121);
imagesc(maxvalnonsym);
colorbar;
title('Max hist val (nonsym)');
subplot(122);
imagesc(maxidxnonsym);
title('Angle idx (nonsym)');

figure(3);
subplot(121);
imagesc(maxval);
colorbar;
title('Max hist val');
subplot(122);
imagesc(maxidx);
title('Angle idx');