%This program demo the analysis and the synthesis of texton features using
%filter banks.
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign

clc;
clear all;
close all;
im_name='SLIM2.bmp';
%Create a synthetic image with 4 classes;
% im=zeros(256,256);
% im(1:128,1:128)=1;
% im(1:128,129:256)=2;
% im(129:256,1:128)=3;
% im(129:256,129:256)=4;

im=im2double(rgb2gray(imread(im_name)));
nrows = size(im,1);
ncols = size(im,2);
%im=im2double((imread(im_name)));

im=im/max(im(:));
im=histeq(im);
%im = imresize(im,0.25);
imagesc(im);%Display the image
title('Original image');
texton_num=20;
texton_calc_en=1;%Enable/Disable texton calculation
[texton,texton_map,filters,f_res,texton_diff]=texton_compute(im,texton_num,'lm',1,'kmean',texton_calc_en);
nscales = 3; ninvscales = 6;
[fe,fo,fs]=computeOrientedEnergyLM(f_res,nscales,ninvscales);
figure
imagesc(texton_map)
title('Texton map');
colorbar
ntexton = size(texton,1);

[texton_image]=compute_texton_pattern(filters,texton);
ntexton1d=ceil(sqrt(ntexton));
figure
imIdx=1;
for rowimageIdx=1:ntexton1d
    for colimageIdx=1:ntexton1d
        if (imIdx<=ntexton)
            subplot(ntexton1d,ntexton1d,imIdx);
            imagesc(texton_image(:,:,imIdx));
            imIdx=imIdx+1;
        end
    end
end

new_f_res=zeros(size(f_res));
for rowimageIdx=1:nrows
    for colimageIdx=1:ncols
        new_f_res(rowimageIdx,colimageIdx,:)=texton(texton_map(rowimageIdx,colimageIdx),:);
    end
end

%Create an indexing image

single_synthetic_map = zeros(nrows,ncols,ntexton);
for textonIdx=1:ntexton
     mask = zeros(size(im));
     idx = find(texton_map==textonIdx);
     mask(idx)=1; %Generate delta function....
     cur_pat = texton_image(:,:,textonIdx);
     single_synthetic_map(:,:,textonIdx)=imfilter(mask,cur_pat,'same');
%     figure
%     imagesc(mask);
%     s=sprintf('Index image %d-th',textonIdx);
%     title(s);
%     colormap gray;
end
figure
%Synthesize the image from its texton map
imagesc(sum(single_synthetic_map,3));
title('Synthetic image from texton map');
colormap gray;

save texton.map texton_image

