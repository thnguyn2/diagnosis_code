%This program is for tesing the idea of single shot learning
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
clc;
clear all;
close all;
addpath(strcat(cd(cd('..')),'\support')); 
datapath = 'H:\QPI data\';
[filenames,glandnames]=findFileName(datapath);
%load the label and the image
lblname = strcat(datapath,'\label\','33C1_resized.mat');
imname = strcat(datapath,'33C1.tif');
textdirhist = strcat(datapath,'\text_dir_hist\','33C1_text_dir_hist.mat');
%Get the directional edgemap

load(textdirhist);
im=imread(imname);
im = imresize(im,size(lblim));
im = im2double(im);
load(lblname);

%Create different patches and corresponding labels
blocksize = 64;
nrows = floor(size(im,1)/blocksize)*blocksize;
ncols = floor(size(im,2)/blocksize)*blocksize;
%Resize both images to match the ass
im = histeq(im(1:nrows,1:ncols));
lblim = lblim(1:nrows,1:ncols);
%Compute the blockidx
[xb,yb]=meshgrid(1:(blocksize):ncols,1:(blocksize):nrows);
nbx = size(xb,1);
nby = size(xb,2);
nblocks = nbx*nby;

blockdata = cell(1,nblocks);
for blockidx=1:nblocks
    blockdata{blockidx} = struct('im',zeros(blocksize,blocksize),'imrot',zeros(blocksize,blocksize),'lbl',zeros(blocksize,blocksize),...
        'maxval',zeros(blocksize,blocksize),'maxidx',zeros(blocksize,blocksize),'mjdir',0,'x',0,'y',0,'mjlbl',0);
end

%Generate a map of regions
regionlbl = zeros(nbx,nby);
mjdirlbl = zeros(nbx,nby);
for bidxx=1:nbx
    for bidxy=1:nby
        blockidx = (bidxy-1)*nbx+bidxx;
        blockdata{blockidx}.im = im(yb(bidxy,bidxx):yb(bidxy,bidxx)+blocksize-1,...
            xb(bidxy,bidxx):xb(bidxy,bidxx)+blocksize-1);
        blockdata{blockidx}.lbl = lblim(yb(bidxy,bidxx):yb(bidxy,bidxx)+blocksize-1,...
            xb(bidxy,bidxx):xb(bidxy,bidxx)+blocksize-1);
        blockdata{blockidx}.x = bidxx;
        blockdata{blockidx}.y = bidxy;
        %Asssign the label of the whole block to the label of the majority
        [tempval,blockdata{blockidx}.mjlbl] = ...
            max(hist(reshape(blockdata{blockidx}.lbl,[1 blocksize^2]),[0 1 2]));
        blockdata{blockidx}.mjlbl = blockdata{blockidx}.mjlbl -1;
        
        blockdata{blockidx}.maxval=maxval(yb(bidxy,bidxx):yb(bidxy,bidxx)+blocksize-1,...
            xb(bidxy,bidxx):xb(bidxy,bidxx)+blocksize-1);
        blockdata{blockidx}.maxidx=maxidx(yb(bidxy,bidxx):yb(bidxy,bidxx)+blocksize-1,...
            xb(bidxy,bidxx):xb(bidxy,bidxx)+blocksize-1);
        %Compute the major direction
        [temp,blockdata{blockidx}.mjdir] = ...
            max(hist(reshape(blockdata{blockidx}.maxidx,[1 blocksize^2]),[1:8]));
        
        regionlbl(bidxx,bidxy)=blockdata{blockidx}.mjlbl;
        mjdirlbl(bidxx,bidxy)=blockdata{blockidx}.mjdir;
        %Compute the direction correspond to the maximum direction
        
    end
end

figure
imagesc(regionlbl);
figure
imagesc(mjdirlbl);


figure(1);
subplot(121);
imagesc(im);
subplot(122);
imagesc(lblim);


