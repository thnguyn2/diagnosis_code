    clc;
    clear all;
    close all;
    addpath(strcat(cd(cd('..')),'\support'));
    labelpath = 'H:\QPI data\label\';
    datapath = 'H:\QPI data\';
    meanstdpath = 'H:\QPI_data\mean&phase';
    [filenames,glandnames]=findFileName(datapath);
    
    %part1: compute the distribution for each image and save it
    h = waitbar(0,'Calculating phase distribution progress...');
    nsamples_per_data = 4000;
    radius_array = [5 10 20 40 80 120 150];
    nradius = length(radius_array);
    npixels = 2048^2;
    nrows = 2048;
    ncols = 2048;
    score_arr_mean = zeros(0,nradius);
    label_arr_mean = zeros(0,nradius);
    score_arr_std = zeros(0,nradius);
    label_arr_std = zeros(0,nradius);

    ys_std = zeros(nradius,1);
    
    colorarr='rbcykmg';    
    sampleidx = 1;
    micronsperpixels = 1/2.8689;
    radius_microns = radius_array*micronsperpixels;
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            label_file_name = strcat(label_name,'.mat');
            disp(['Calculating phase histogram: ' label_name ' ...']);
            org_im = imread(cur_file_name); %Load the original data
            load(strcat(labelpath,label_file_name));
            org_im = imresize(org_im,[2048 2048],'nearest');
            lblim = imresize(lblim,[2048 2048],'nearest');
            org_im = cast(org_im,'single');
            glandidxorg = find(lblim==1);
            stromaidxorg = find(lblim==2);
            
            %Compute ROC curve for other radius
            pixelidx=[1:npixels];
            [y_coord,x_coord]=ind2sub([2048 2048], pixelidx); %Row and column indices of current pixels
            y_coord = y_coord(:);
            x_coord = x_coord(:);
            intphaseIm = integralImage(org_im);
            intphaseIm=intphaseIm(1:end-1,1:end-1);           
            
            intphaseSqIm=integralImage(org_im.^2);
            intphaseSqIm = intphaseSqIm(1:end-1,1:end-1);
            stdim = zeros(nrows,ncols);
          
           
            for radidx = 1:nradius
                radius  = radius_array(radidx);
                x_offset = [-(radius/2) (radius/2) (radius/2) -(radius/2)];
                y_offset = [-(radius/2) -(radius/2) (radius/2) (radius/2)];
                x_offset = floor(x_offset);
                y_offset = floor(y_offset);
                neigh_x = repmat(x_coord,[1 4])+ones(npixels,1)*x_offset;
                neigh_y = repmat(y_coord,[1 4])+ones(npixels,1)*y_offset;
                %Check for boundary conditions
                idx = find(neigh_x<1);
                %Check for boundary conditions
                idx = find(neigh_x<1);
                neigh_x(idx)=1;
                idx = find(neigh_x>ncols);
                neigh_x(idx)=ncols;
                idx = find(neigh_y<1);
                neigh_y(idx)=1;
                idx = find(neigh_y>nrows);
                neigh_y(idx)=nrows;
                %Convert into linear indexing....
                neigh_coord=sub2ind([nrows ncols],neigh_y,neigh_x);
                %Compute the number of neighbors for each pixels
                nneighbors = (neigh_x(:,2)-neigh_x(:,1)).*(neigh_y(:,3)-neigh_y(:,2));
                averagephase = zeros(nrows,ncols);
                averagephaseSum = zeros(nrows,ncols);
                averageSqphase = zeros(nrows,ncols);
                %Compute the average phase in a window
                
                averagephaseSum(pixelidx) =  intphaseIm(neigh_coord(:,3))+ intphaseIm(neigh_coord(:,1))-...
                         intphaseIm(neigh_coord(:,4))- intphaseIm(neigh_coord(:,2));
                %Compute the standard deviation image
                averageSqphase(pixelidx) =  intphaseSqIm(neigh_coord(:,3))+ intphaseSqIm(neigh_coord(:,1))-...
                         intphaseSqIm(neigh_coord(:,4))- intphaseSqIm(neigh_coord(:,2));
                %Compute the mean phase image 
                averagephase(:) = averagephaseSum(:)./nneighbors;                
                stdim = ((averageSqphase - 2*averagephaseSum.*averagephase)./reshape(nneighbors,[2048 2048])+averagephase.^2);
                stdim(:) = sqrt(stdim(:));   
                
                %Add the mean phase for evaluation
                glandphase = averagephase(glandidxorg); %Compute the score for gland and stroma
                stromaphase = averagephase(stromaidxorg);
                glandphase = glandphase(1:2*floor(size(glandphase,1)/nsamples_per_data):end);
                stromaphase = stromaphase(1:2*floor(size(stromaphase,1)/nsamples_per_data):end);
                score = [glandphase(:);stromaphase(:)];
                labeldata=[ones(size(glandphase));2*ones(size(stromaphase))];
                score_arr_mean(sampleidx:sampleidx + length(score)-1,radidx)=score;
                label_arr_mean(sampleidx:sampleidx + length(score)-1,radidx)=labeldata;
                
                %Add the std of the phase for evaluation
                glandstd = stdim(glandidxorg); %Compute the score for gland and stroma
                stromastd = stdim(stromaidxorg);
                glandstd = glandstd(1:2*floor(size(glandstd,1)/nsamples_per_data):end);
                stromastd = stromastd(1:2*floor(size(stromastd,1)/nsamples_per_data):end);
                score = [glandstd(:);stromastd(:)];
                score_arr_std(sampleidx:sampleidx + length(score)-1,radidx)=score;
                label_arr_std(sampleidx:sampleidx + length(score)-1,radidx)=labeldata;
                
                
                
                [xs_std,ys_std,t,auc_std]=perfcurve(label_arr_std(:,radidx),score_arr_std(:,radidx),2);
                figure(1);
                plot(xs_std,ys_std,colorarr(radidx),'linewidth',2);
                hold on
                figure(2);
                [xs_mean,ys_mean,t,auc_mean]=perfcurve(label_arr_mean(:,radidx),score_arr_mean(:,radidx),2);
                figure(2);
                plot(xs_mean,ys_mean,colorarr(radidx),'linewidth',2);
                hold on;
                
                figure(3);
                subplot(121);
                averagephase = averagephase/65536*4.0-0.5;
                imagesc(averagephase);
                colorbar
                title('Mean phase');
                subplot(122);
                stdim = stdim/65536*4.0-0.5;
                imagesc(stdim);
                colorbar;
                title('Phase std');
            end
            sampleidx = sampleidx + length(score);
            figure(1);
            hold off;
            legend('1.74 um','3.49 um','6.97 um','13.94 um','27.88 um','41.83 um','52.28 um');
            title('ROC curve for the mean phase');
            figure(2);
            hold off;
            legend('1.74 um','3.49 um','6.97 um','13.94 um','27.88 um','41.83 um','52.28 um');
            title('ROC curve for the phase std');
        end
    end
   save('avgstd.mat','score_arr_mean','score_arr_std','label_arr_mean','xs_mean','ys_mean','xs_std','ys_std');
    
%Draw the ROC curve for mean and std of the phase. This section of the code
%is for saving the ROC curve
load('avgstd.mat');
nradius = 7;
figure(1);
ys_std = cell(nradius,1);
xs_std = cell(nradius,1);
ys_mean = cell(nradius,1);
xs_mean = cell(nradius,1);
for radidx = 1:nradius
    [xs_std{radidx},ys_std{radidx},t,auc_std]=perfcurve(label_arr_mean(:,radidx),score_arr_std(:,radidx),2);
    figure(1);
    plot(xs_std{radidx},ys_std{radidx},'b','linewidth',2);
    hold on
    [xs_mean{radidx},ys_mean{radidx},t,auc_mean]=perfcurve(label_arr_mean(:,radidx),score_arr_mean(:,radidx),2);
    figure(2);
    plot(xs_mean{radidx},ys_mean{radidx},'g','linewidth',2);
end
save('avgstd.mat','xs_mean','ys_mean','xs_std','ys_std','-append');
close(h);
