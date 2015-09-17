clc;
clear all;
close all;
load texton.mat;
im_name='SLIM4.bmp';

im=im2double(rgb2gray(imread(im_name)));
nrows = size(im,1);
ncols = size(im,2);
ntexton=size(texton_image,3);
%Compute the correlation between the input image with the texton pattern
f_cor = zeros(nrows,ncols,ntexton);
disp('Computing correlation with texton image...');
for textonIdx=1:ntexton
    cur_pat = texton_image(:,:,textonIdx);
    f_cor(:,:,textonIdx)=imfilter(im,cur_pat,'same');
end
[max_val,texton_map]=max((f_cor),[],3);
figure
imagesc(im);
colormap gray;
figure
imagesc(texton_map);
colormap jet
colorbar
single_synthetic_map = zeros(nrows,ncols,ntexton);
for textonIdx=1:ntexton
     mask = zeros(size(im));
     idx = find(texton_map==textonIdx);
     mask(idx)=1; %Generate delta function....
     cur_pat = texton_image(:,:,textonIdx);
     single_synthetic_map(:,:,textonIdx)=imfilter(mask,cur_pat,'same');
end
imagesc(sum(single_synthetic_map,3));
title('Synthetic image from texton map');