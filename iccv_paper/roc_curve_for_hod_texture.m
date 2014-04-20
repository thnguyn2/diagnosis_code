%This function evaluate the ROC curve of the HOD feature
clc;
clear all;
close all;
%This code display different ROC curves for different features
%For SVM on textons
[filenames,glandnames]=findFileName('H:\QPI data\');
histpath = 'H:\QPI data\text_dir_hist\';
h = waitbar(0,'Reading data textons...');
nfilestotal = 0;
f_total = zeros(0,1);
gt_total = zeros(0,1); %Groundtruth
color_arr='rbgm';



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
           texture_dir_file_name = strcat(label_name,'_lm_fr.mat'); %For loading the texture direction
           texture_dir_hist_file_name =  strcat(label_name,'_text_dir_hist.mat');
           texton_idx_file_name = strcat(label_name,'_texton.mat'); %For loading the label
           
           disp(['Adding SVM data of: ' label_name ' ...']);
           load(strcat(histpath,texture_dir_hist_file_name),'klsym','klnonsym','lblim');
           lblim = imresize(lblim,[256 256],'nearest');
           klsym = imresize(klsym,[256 256],'nearest');
           klnonym = imresize(klnonsym,[256 256],'nearest');
           
           glandidx = find(lblim==1);
           stromaidx = find(lblim==2);
           klgland = klsym(glandidx);
           klstroma = klsym(stromaidx);
           label_arr=[ones(length(glandidx(:)),1);(-1)*ones(length(stromaidx(:)),1)];
           value_arr = [klgland(:);klstroma(:)];
           [xssym,yssym]=perfcurve(label_arr,value_arr,-1);
           figure(1);
           plot(xssym,yssym,color_arr(classidx));
           hold on;
           xlabel('False positive rate (stroma)[sym]'); ylabel('True positive rate (gland)[sym]');
           
           
           klgland = klnonsym(glandidx);
           klstroma = klnonsym(stromaidx);
           label_arr=[ones(length(glandidx(:)),1);(-1)*ones(length(stromaidx(:)),1)];
           value_arr = [klgland(:);klstroma(:)];
           [xsnon,ysnon]=perfcurve(label_arr,value_arr,-1);


           figure(2);
           plot(xsnon,ysnon,color_arr(classidx));
           hold on;
           xlabel('False positive rate (stroma)[nonsym]'); ylabel('True positive rate (gland)[nonsym]');

     end
end



