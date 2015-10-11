%This is the demonstration for using the classification code
%To segment the image, put it in the val_data as below, make sure the pixel
%ratio is correct

%Testing part;
clc;
clear all;
close all;
imname = 'H:\Cores and ROIs for Tans trainer [Hassaan]\New folder\val_data\A0.tif'; %Make sure the image is 16 bit with some background in the topleft corner
nscales = 4; ninvscales = 5; ndir=6;
windowsize = 60;
nrows = 2048; ncols = 2048; %Dimensions after downsampling
phaseim = single(imresize(imread(imname),[nrows ncols]));
figure(1);
subplot(221);imagesc(phaseim);colormap jet;colorbar; title('Phase image to be segmented');

%First steps - Compute filter response
[mr_data,filters,L]=compute_mr_response(phaseim,nscales,ninvscales,ndir);

%Load the extracted texton data
load(strcat(pwd,'\texdir\','kmeans_res_50_clusters.mat'));

%Second stage: compute the vector quantization on
new_text_map = vect_quant(mr_data,texton,nrows,ncols,0); %Vector quantization on the texton
subplot(222);imagesc(new_text_map);colormap jet;title('Texton indexing map');

%Third stage: compute the histogram of textons
ntextons = size(texton,1);
histim = computehistogramoftexton(new_text_map,ntextons,windowsize);

%Forth stage: specify the no-lumen region and run the classification on
%these area
ds_mean=mean(mean(phaseim(1:30,1:30))); %Find the phase of the background region
mask = (phaseim>1.08*ds_mean); %Make sure what we have is good for training
se = strel('disk',5);
newmask = imclose(mask,se);
newinvmask = imclose(bwareaopen(imcomplement(newmask),5000),se); %Guess what i am doing :)
newmask = bwareaopen(imclose(imcomplement(newinvmask),se),5000);
nonlumenidx= find(newmask==1); %This is the coordinates of the pixels we are going to classify

%Step 5: run the classifier on the non-lumen pixels
labelim=zeros(nrows,ncols);
load(strcat(pwd,'\Classifier\rfstruct_ws60_May_8th_2015_2_classes.mat'));%load the random forest classifier

%Convert 3D dataset into 2D dataset, each column is one band. Each row is one pixel
testdata = zeros(nrows*ncols,ntextons);
for bandidx=1:ntextons
    curband = histim(:,:,bandidx);
    curvect = curband(:);
    testdata(:,bandidx)=curvect(:);
end
valset = testdata(nonlumenidx,:);
nevalsample = size(valset,1);
nsampleperbatch =80000;
nbatch=ceil(size(valset,1)/nsampleperbatch);
outputlabel=zeros(size(nonlumenidx));
for batchidx=1:nbatch
    disp(['Batch idx ', num2str(batchidx) '/' num2str(nbatch)]);
    cures=rfstruct.predict(valset((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,nevalsample),:));
    outputval((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,nevalsample),:)=cures;
end
putputval = str2num(cell2mat(outputval);
                      
    


