%This source code perform multi-resolution segmentation for glands vs
%stroma. First, it will swipe through all folder, look for two fuzzy map of
%class prob at two different radii. Combine it and produce a final
%segmentation map
clc;
clear all;
close all;
datapath = 'H:\TMA_cores_and_diagnosis\';
texton_dir=strcat(datapath,'texdir\');
label_dir=strcat(datapath,'label\2class_results\');
addpath(strcat(cd(cd('..')),'\support'));
[filenames,glandnames,classlist,alltiffilelist]=findFileNameFromROIs(datapath);
radius = 90;
retrain = 1;
param.glandthresh = 0.5; %Threshold for gland classification
param.boundthreshforrejoin = 0.2;
param.mindistancebetweenglands = 15;
param.glandboundmin = 20;
param.minlumenarea =1500;
param.minglandsize = 5000;
for fileidx = 1:length(alltiffilelist)
                    cur_file_name = alltiffilelist{fileidx};
                    dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                    slash_pos = strfind(cur_file_name,'\');
                    label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                    disp(['Working on ' label_name ', Sample Idx: ' num2str(fileidx)])
                    seg_im_name_prob = strcat(label_dir,label_name,'_',num2str(radius),'_2cl.tif');  %This is the probability map
                    seg_im_name = strcat(label_dir,label_name,'_seg_multi_res.tif'); %Name of the class map
                    if ((exist(seg_im_name_prob)&(~exist(seg_im_name)))|(retrain==1))
                        fmap=imread(seg_im_name_prob); %Read the first map
                        lblim=segmentation_refine(fmap,param);
                        writeTIFF(lblim,seg_im_name);
                    end               
end        
                 
