%This program demo the analysis and the synthesis of texton features using
%filter banks.
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign

clc;
clear all;
close all;
im_name='SLIM6.bmp';
%Create a synthetic image with 4 classes;
im=zeros(256,256);
im(1:128,1:128)=1;
im(1:128,129:256)=2;
im(129:256,1:128)=1;
im(129:256,129:256)=2;

im=im2double(rgb2gray(imread(im_name)));
%im=im2double((imread(im_name)));
im=im/max(im(:));
%im=histeq(im);
%im = imresize(im,0.25);
imagesc(im);%Display the image
title('Original image');
colormap gray;

[f_res,mage,ee,eo,e,dir_map]=filter_response_compute_with_conf(im);
figure
subplot(131)
imagesc(im)
subplot(132)
imagesc(sum(eo,3))
subplot(133)
imagesc(sum(ee,3))

figure
imagesc(sum(eo,3)+sum(ee,3));
title('Combined boundary map')
figure
imagesc(sum(eo,3));
title('Odd boundary map');
figure
imagesc(sum(ee,3));
title('Even boundary map');


imedge=edge(im,'canny');%Find the edge image
