%This code overlays the grayscale image with the RGB image. The grayscale
%image will serve as the modulation to the RGB channel
clear all;
close all;
clc;
labelfilename = 'C:\Users\thnguyn2\OneDrive\optics_materials\Cancer diagnosis\Figures used for papers\O2_seg_multi_res_3072_strict_hf.tif';
grayscalefilename = 'C:\Users\thnguyn2\OneDrive\optics_materials\Cancer diagnosis\Figures used for papers\O2_small.tif';
labelim = imread(labelfilename);
grayim = imread(grayscalefilename);
grayim = cast(grayim/256,'single');
b = cast((labelim==0)*255,'single'); %Lumen
g = cast((labelim==1)*255,'single'); %Gland
r = cast((labelim==2)*255,'single'); %Stroma
[nrows,ncols]=size(grayim);
overlaidim = zeros(nrows,ncols,3);
overlaidim(:,:,1)=cast(r.*(grayim)/256^2,'single');
overlaidim(:,:,2)=cast(g.*(grayim)/256^2,'single');
overlaidim(:,:,3)=cast(b.*(grayim)/256^2,'single');
overlaidim = overlaidim*5; %Make the dynamic range better
idx = find(overlaidim>1.0);
overlaidim(idx)=1.0;
finalim = 0.21*overlaidim(:,:,1) + 0.71*overlaidim(:,:,2) + 0.07*overlaidim(:,:,3);
figure(1);
imagesc(labelim);
figure(2);
imagesc(grayim);
figure(3);
imagesc(overlaidim);truesize;