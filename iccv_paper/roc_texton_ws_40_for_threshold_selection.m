clc;
clear all;
close all;
%This code display different ROC curves for different samples. This will
%help us determine the best threshold
%For SVM on textons
[filenames,glandnames]=findFileName('H:\QPI data\');
svm_dir = 'H:\QPI data\svm\texton_ws40_512\';
h = waitbar(0,'Reading data textons...');
nfilestotal = 0;
rocnpoints = 679134;
f_total = zeros(0,1);
gt_total = zeros(0,1); %Groundtruth
color_arr='rbgm';
load('svmstruct.mat');



for classidx=1:4 %Go through different classes
    nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
    for sampleIdx=1:nsamples
           nfilestotal =nfilestotal+1; 
           waitbar(sampleIdx/nsamples,h,'Progress...')
           cur_file_name = filenames{classidx,1}{sampleIdx,1};
           %Check to see if the label file exist
           dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
           slash_pos = strfind(cur_file_name,'\');
           label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
           svm_file_name = strcat(label_name,'_svm_200k_ws40.mat');          
           disp(['Adding SVM data of: ' label_name ' ...']);
           load(strcat(svm_dir,svm_file_name));
           
           %Compute the ROC curve, threshold, auc
           glandidx = find(lblim==1);
           stromaidx = find(lblim==-1);
           ngland = length(glandidx(:));
           nstroma = length(stromaidx(:));
           fimgland = fim(glandidx);
           fimstroma =fim(stromaidx);
           score = [fimgland(:);fimstroma(:)];
           labelmap = [ones(ngland,1);2*ones(nstroma,1)];
           [xs,ys,t,auc] = perfcurve(labelmap,score,1);
           
           figure(1);
           plot(xs,ys,color_arr(classidx));
           hold on;
           xlabel('P_F (stroma)'); ylabel('P_D (gland)');
         
           figure(5);
           plot(t,xs,color_arr(classidx),'linewidth',2);
           xlabel('Threshold');
           ylabel('P_F (stroma)');
           hold on;
           
           figure(6);
           plot(t,1-ys,color_arr(classidx),'linewidth',2);
           xlabel('Threshold');
           ylabel('P_F (gland)');
           hold on;
           
           figure(2);
           imagesc(fim);
           title('Functional distance');
           
           newlblim = zeros(size(lblim));
           thresh = -0.028;
           candidateidx = find(lblim~=0);
           candidatelabel = (-1)*ones(size(candidateidx));
           glandidx = find(fim(candidateidx)>thresh);
           candidatelabel(glandidx)=1;
           newlblim(candidateidx)=candidatelabel;
           lumenidx = find(newlblim==0);
           newlblim(lumenidx)=-1;
           %Apply image dilation
           se = strel('disk',1);
           idxim = im2bw(newlblim,0);
           idxim = imclose(idxim,se);
           newlblim = cast(idxim,'single')*2.0-1.0;
           newlblim(lumenidx)=0;
           
           colorbar
           figure(3);
           imagesc(newlblim);
           title('Labeling results from texton')
           figure(4);
           imagesc(lblim);
           title('Human labeling');
         %  a=1;
    end
end



