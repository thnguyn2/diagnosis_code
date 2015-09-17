clc;
clear all;
close all;
dataname = 'H:\QPI data\ForNuclei\33A8_lm_fr_fine_scale.mat';
imname = 'H:\QPI data\33A8.tif';
load(dataname,'edir','emap','emag','f_res');
edgemap = emap(:,:,1);
dirmap = edir(:,:,1);
ncols = 2048;
nrows = 2048;
ndirs = 20; %Number of pixel direction


mask = im2bw(edgemap,0.5);
s = regionprops(mask,'Area','PixelList');
area_arrays = cat(1,s.Area);
rejectidx = find(area_arrays<=3);%Get rid of the clutter
newedgemap = edgemap;
newdirmap = dirmap;
for idx = 1:length(rejectidx)
   curcoord=s(rejectidx(idx),1).PixelList;
   coord1d = sub2ind([nrows ncols],curcoord(:,2),curcoord(:,1));
   newedgemap(coord1d)=0;
   newdirmap(coord1d)=0;
end
edgemap = newedgemap;
dirmap = newdirmap;
% pixelcoord = cat(1,s.Centroid);
% centroids = round(centroids);
% centroids1d = sub2ind([nrows ncols],centroids(:,2),centroids(:,1));



% edgemap = zeros(2048,2048);
% dirmap = zeros(2048,2048);
% edgemap(end/2,end/2)=1;
% dirmap(end/2,end/2)=3;


diredgemap = edgemap.*dirmap;
diredgemap = (diredgemap-1)/20*180; %(Angle here is in degree)
idx = find(diredgemap<0);
diredgemap(idx)=0;

figure(1);
imagesc(diredgemap);
colorbar
title('Edge direction map');
figure(6);
imagesc(diredgemap(end/4:3*end/4,end/4:3*end/4));
colorbar;
title('Dir edgemap zoomed');

%Computing Hough map for each direction
edgeidx = find(edgemap==1);
%Define a radius where all pixels in this radius will be 
minrad = 8;
maxrad = 14;
radval = minrad:maxrad;
radval = radval(:)';
npixels = length(edgeidx);
radval = repmat(radval,[npixels 2]);

[y_coord,x_coord]=ind2sub([2048 2048], edgeidx); %Row and column indices of the edge pixels
angleval = 180*(dirmap(edgeidx)-1)/ndirs;
angleval = angleval(:);
angleval = repmat(angleval,[1 size(radval,2)/2]); %Generate the first half of the angle
angleval = [angleval (angleval+180)];
%Now, compute the direction of the points that we should add

x_offset = -sind(angleval).*radval;
y_offset = cosd(angleval).*radval;
nnbrs = size(x_offset,2);
neigh_x =repmat(x_coord,[1 nnbrs]) + x_offset;
neigh_y =repmat(y_coord,[1 nnbrs]) + y_offset;
%Check for boundary conditions
idx = find(neigh_x<1);
neigh_x(idx)=1;
idx = find(neigh_x>ncols);
neigh_x(idx)=ncols;
idx = find(neigh_y<1);
neigh_y(idx)=1;
idx = find(neigh_y>nrows);
neigh_y(idx)=nrows;
neigh_x = round(neigh_x);
neigh_y = round(neigh_y);

disp('Computing Hough map data...');
houghmap = zeros(nrows,ncols);
dirhoughmap = zeros(nrows,ncols,ndirs); %This is the directional Hough map
%Now, go for each neighbor and add up the count
blob_std =5;
hlog=-fspecial('log',[30 30],blob_std);

for nbidx = 1:nnbrs
    nbidx
    votemask = zeros(nrows,ncols);
    coord = [neigh_x(:,nbidx) neigh_y(:,nbidx)];
    [coordval] = unique(coord,'rows');%Find the unique row elements
    [temp1,temp2,idx] = unique([coordval;coord],'rows');%Find the unique row elements
    npixelsassigned = size(coordval,1);
    %Now, truncate the first part of the indexing so that only the
    %overlapping's histogram is calculated
    idx = idx(npixelsassigned+1:end);
    nvotes = histc(idx,1:npixelsassigned);
    updatedpixelidx = sub2ind([nrows ncols],coordval(:,2),coordval(:,1));
    votemask(updatedpixelidx)=nvotes;
    %Suppress falty information at the edge
    votemask(1,:)=0;
    votemask(:,1)=0;
    votemask(end,:)=0;
    votemask(:,end)=0;
    houghmap = houghmap + votemask;
    figure(5);
    %imagesc(votemask(3*end/8:5*end/8,3*end/8:5*end/8));
    imagesc(votemask);
    colorbar
    title('Vote map');
    figure(2);
    %imagesc(houghmap(3*end/8:5*end/8,3*end/8:5*end/8));
    imagesc(houghmap);
    colorbar;
    title('Hough map');
     drawnow;
   
end
con_hough_map=imfilter(houghmap,hlog,'same');

figure(2);
imagesc(con_hough_map);
colorbar;
title('Concentration hough map');
%Find the index of the center
idx = find(con_hough_map>0.2*max(con_hough_map(:)));
 
%Now, suppress the regions of blobs and replace it with the centrold
mask = zeros(nrows,ncols);
mask(idx)=1;

%Before moving on, elminiate all regions with two few pixels
mask = im2bw(mask,0.5);
s = regionprops(mask,'Area','PixelList');
area_arrays = cat(1,s.Area);
rejectidx = find(area_arrays<=5);%Get rid of two small center
for idx = 1:length(rejectidx)
   curcoord=s(rejectidx(idx),1).PixelList;
   coord1d = sub2ind([nrows ncols],curcoord(:,2),curcoord(:,1));
   mask(coord1d)=0;
end
oldmask = mask;
mask = im2bw(mask,0.5);
s = regionprops(mask,'centroid');
centroids = cat(1,s.Centroid);
centroids = round(centroids);
centroids1d = sub2ind([nrows ncols],centroids(:,2),centroids(:,1));
%centroids1d = idx;
edgemap(centroids1d)=2;

figure(3);
imagesc(edgemap);
colorbar;
title('Edge map');
%Next, call the phase image and compute the 
im = imread(imname);
im = cast(im,'single');
im = imresize(im,[2048 2048]);
im(centroids1d)=80000;
figure(4)
imagesc(im);
colorbar
colormap gray

