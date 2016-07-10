%This source code calculate the Fractal dimension
%Author: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
function demo
clc;
clear all;
close all;
I = imread('http://4.bp.blogspot.com/-aHCfmDvyzFU/Un_U-Neo_GI/AAAAAAAAGpQ/DWzjztkh4HM/s1600/sierpinski.png');
Ibw = ~im2bw(I); %Note that the background need to be 0
figure(1);
imagesc(Ibw);
colormap gray;
tic
dim_val=BoxCountfracDim(Ibw) %Compute the box-count dimension
toc
tic
D = hausDim(Ibw)
toc
%Convert the image toe 
end

function dim_val=BoxCountfracDim(im)
%Calculate the fractal dimension
%Input:
%   im: input image
%Output:
%   dim_val: value of the fractal dimension
%Note:
%    Returns the Haussdorf fractal dimension D of an object represented by the
%    binary image I. Nonzero pixels belong to an object and 0 pixels 
%    constitute the background.
%
%    Algorithm
%    ---------
%    1 - Pad the image with background pixels so that its dimensions are a 
%        power of 2.
%    2 - Set the box size 'e' to the size of the image.
%    3 - Compute N(e), which corresponds to the number of boxes of size 'e' 
%        which contains at least one object pixel.
%    4 - If e > 1 then e = e / 2 and repeat step 3.
%    5 - Compute the points log(N(e)) x log(1/e) and use the least squares 
%        method to fit a line to the points.
%    6 - The returned Haussdorf fractal dimension D is the slope of the line.
%Author: Tan H. Nguyen
    %First, padd the image with 0 to make its dimension a power of 2
    newdim = 2^max(nextpow2(size(im)));
    im_pad = padarray(im,[newdim-size(im,1) newdim-size(im,2)],0,'post');
    figure(1);
    imagesc(im_pad);
    maxpow2 = log2(newdim);
    boxCount_arr = zeros(1,maxpow2+1);
    res_arr = zeros(1,maxpow2+1);
    curboxsize = newdim;
    boxCount_arr(1)=1;
    res_arr(1)=1.0/curboxsize;
    curim = im_pad;
    for idx = 1:maxpow2
        sum_val = sum(curim(:));%Go with the finest scale to the coarsest
        boxCount_arr(2+maxpow2-idx)=sum_val;
        res_arr(2+maxpow2-idx)=2^(-(idx-1));
        curim = ((curim(1:2:end,1:2:end)+curim(1:2:end,2:2:end)+curim(2:2:end,1:2:end)+curim(2:2:end,2:2:end))>0);
    end
    D = polyfit(log(res_arr), log(boxCount_arr), 1);
    dim_val = D(1);
    
end


