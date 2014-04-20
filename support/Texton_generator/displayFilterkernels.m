clc;
clear all;
close all;
F = makeLMfilters;
a1 = F(:,:,17);
a2 = F(:,:,18);
a3 = F(:,:,19);
a4 = F(:,:,20);
a5 = F(:,:,21);
a6 = F(:,:,22);
a7 = F(:,:,23);
a8 = F(:,:,24);

b1 = F(:,:,57);
b2 = F(:,:,58);
b3 = F(:,:,59);
b4 = F(:,:,60);
b5 = F(:,:,61);
b6 = F(:,:,62);
b7 = F(:,:,63);
b8 = F(:,:,64);

c1 = F(:,:,81);
c1 = c1/max(c1(:));
c2 = F(:,:,83);
c2 = c2/max(c2(:));
c3 = F(:,:,85);
c3 = c3/max(c3(:));
c4 = F(:,:,87);
c4 = c4/max(c4(:));
c5 = F(:,:,82);
c5 = c5/max(c5(:));
c6 = F(:,:,84);
c6 = c6/max(c6(:));
c7 = F(:,:,86);
c7 = c7/max(c7(:));
c8 = F(:,:,88);
c8 = c8/max(c8(:));

figure(1);
imagesc([   a1 a2 a3 a4 a5 a6 a7 a8;...
            b1 b2 b3 b4 b5 b6 b7 b8]);
colormap jet;
axis off;
truesize;
figure(2);
subplot(211);
imshow([c1 c2 c3 c4]);
subplot(212);
imagesc([c5 c6 c7 c8]);
colormap jet;
truesize;
axis off;

%Display the texton images---
load('texton_data.mat');
dir = 1;
%indexset = [0+dir:8:80];
indexset = [81:90];
%indexset = [indexset(:);[81:90]'];
% newfilterset = zeros(size(F,1),size(F,2),length(indexset));
% newfilterset(:,:,:) = F(:,:,indexset);
% texton_image = compute_texton_pattern(newfilterset,fin_texton);

newfilterset = zeros(size(F,1),size(F,2),length(indexset));
newfilterset(:,:,:) = F(:,:,indexset);
texton_image = compute_texton_pattern(newfilterset,fin_texton(:,1:length(indexset)));

%display the imagesc
nimages = size(texton_image,3);
figure(1);
nsqr= ceil(sqrt(nimages));
for imidx =1:nimages
    subplot(7,8,imidx);
    imagesc(texton_image(:,:,imidx));
    colormap gray;
    axis off;
 end
figure(2);
for textonidx=1:size(fin_texton,1)
    plot(fin_texton(textonidx,:));
    hold on;
end