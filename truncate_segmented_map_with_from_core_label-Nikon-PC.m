function truncate_segmented_map_with_from_core_label
    clc;
    clear all;
    close all;
    datapath = 'H:\TMA_cores_and_diagnosis\';
    texton_dir=strcat(datapath,'texdir\');
    label_dir=strcat(datapath,'label\2class_results\');
    core_groundtruth_dir = strcat(datapath,'diagnosis_of_vicky\');%This is the groundtruth for the core
    addpath('.\support')
    addpath(strcat(cd(cd('..')),'\support'));
    [filenames,glandnames,classlist,alltiffilelist]=findFileNameFromROIs(datapath);
    savephaseim = 0;
    phaseim = -1;
    for fileidx = 1:length(alltiffilelist)
        cur_file_name = alltiffilelist{fileidx};
        %cur_file_name = 'E:\TMA_cores_and_diagnosis\d3+3\J4.tif';
        dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
        slash_pos = strfind(cur_file_name,'\');
        label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
        tic;
            disp(['Working on ' label_name ', Sample Idx: ' num2str(fileidx)])
            if (savephaseim)
                phaseim = imread(strcat(cur_file_name(1:end-4),'.tif'));%Downsampled image  
            end
            labelim = imread(strcat(label_dir,label_name,'_seg_multi_res.tif')); %2 class map.
            %Load the core diagnosis and use it to crop the segmented map
            coregt_map_file_name = strcat(core_groundtruth_dir,label_name,'_roi_diag.tif');
            roi_loc_file_name = strcat(core_groundtruth_dir,label_name,'_roi_diag_loc.mat');
            
            core_index_texton_name = strcat(texton_dir,label_name,'_texton_index_map.tif');
            if (exist(coregt_map_file_name))
                load(roi_loc_file_name);%Load the position and the class for all the rois
                texidxmap = imread(core_index_texton_name);%Load the histogram of texton
                create_crop_imagesc(3,roiindexlist,roitypelist,phaseim,texidxmap,labelim,label_name,datapath);%Create all the image of grade 3 cores
                create_crop_imagesc(4,roiindexlist,roitypelist,phaseim,texidxmap,labelim,label_name,datapath);%Grade 4 case
                create_crop_imagesc(5,roiindexlist,roitypelist,phaseim,texidxmap,labelim,label_name,datapath);%Grade 5 case
                create_crop_imagesc(-1,roiindexlist,roitypelist,phaseim,texidxmap,labelim,label_name,datapath);%Normal case
                create_crop_imagesc(-3,roiindexlist,roitypelist,phaseim,texidxmap,labelim,label_name,datapath);%BPH case
                 create_crop_imagesc(-5,roiindexlist,roitypelist,phaseim,texidxmap,labelim,label_name,datapath);%Grade 5 case
               

            else
                disp(['Core ' label_name ' ground truth does not exist...'])
            end
            tproc = toc;
            disp(['Processing time: ' num2str(tproc) ' (s)']);

      end
end

function [tempphaseim,temptextonidx]=create_crop_imagesc(class,roiindexlist,roitypelist,phaseim,texidxmap,labelim,label_name,datapath)
     g3folder = strcat(datapath,'diagnosis_of_vicky\g3\');
     g4folder = strcat(datapath,'diagnosis_of_vicky\g4\');
     g5folder = strcat(datapath,'diagnosis_of_vicky\g5\');
     nmfolder = strcat(datapath,'diagnosis_of_vicky\nm\');
     bphfolder = strcat(datapath,'diagnosis_of_vicky\bph\');
     hgpfolder = strcat(datapath,'diagnosis_of_vicky\hgp\');
     if (phaseim~=-1)
        dimratio  = size(phaseim,1)/size(labelim,1);
     end
     foundcoreidx = find(roitypelist==class); %Find all the rois of the current class
     if (~isempty(foundcoreidx))
         nregions = length(foundcoreidx); 
         for regionidx=1:nregions
               curidx1 = roiindexlist{foundcoreidx(regionidx)}; 
               curarea = length(curidx1);
               if (curarea>=1000)
                   [currowidx,curcolidx]=ind2sub(size(labelim),curidx1);
                   r1 = min(currowidx);
                   c1 = min(curcolidx);
                   r2 = max(currowidx);
                   c2 = max(curcolidx);
                   bbnrows = r2-r1+1;
                   bbncols = c2-c1+1;
                   labelmap = zeros(bbnrows,bbncols);
                   temptextonidx = zeros(bbnrows,bbncols);
                   currowidx = currowidx+1-r1;
                   curcolidx = curcolidx+1-c1;
                   curidx = currowidx + (curcolidx-1)*bbnrows;
                   labelmap(curidx) = labelim(curidx1);
                   temptextonidx(curidx) = texidxmap(curidx1);
                  
                   if (phaseim~=-1) %If we want to save the phase image as well
                       tempphaseim = phaseim(round(r1*dimratio):round(r2*dimratio),round(c1*dimratio):round(c2*dimratio));
                       tempphaseim = tempphaseim.*cast((imresize(temptextonidx,size(tempphaseim),'nearest')~=0),'uint16');
                   end
                   if (class==4)
                         %writeTIFF(tempphaseim,strcat(g4folder,label_name,'_small_',num2str(regionidx),'_g4.tif'));%Save the phase image
                         %writeTIFF(tempphaseim,strcat(g4folder,label_name,'_',num2str(regionidx),'_g4.tif'),'uint16');%Save the phase image
                         %writeTIFF(labelmap,strcat(g4folder,label_name,'_lbl_',num2str(regionidx),'_g4.tif'));%Save the phase image
                         %writeTIFF(temptextonidx,strcat(g4folder,label_name,'_texidx_',num2str(regionidx),'_g4.tif'));%Save the phase image
                         save(strcat(g4folder,label_name,'_coord_',num2str(regionidx),'_g4.mat'),'r1','r2','c1','c2');
                   elseif (class==3)
                         %writeTIFF(tempphaseim,strcat(g3folder,label_name,'_small_',num2str(regionidx),'_g3.tif'),'uint16');%Save the phase image
                         
                         %writeTIFF(tempphaseim,strcat(g3folder,label_name,'_',num2str(regionidx),'_g3.tif'),'uint16');%Save the phase image
                         %writeTIFF(labelmap,strcat(g3folder,label_name,'_lbl_',num2str(regionidx),'_g3.tif'));%Save the phase image
                         %writeTIFF(temptextonidx,strcat(g3folder,label_name,'_texidx_',num2str(regionidx),'_g3.tif'));%Save the phase image
                         save(strcat(g3folder,label_name,'_coord_',num2str(regionidx),'_g3.mat'),'r1','r2','c1','c2');
           
                   elseif (class==5)
                       % writeTIFF(tempphaseim,strcat(g5folder,label_name,'_small_',num2str(regionidx),'_g5.tif'),'uint16');%Save the phase image
                         
                       % writeTIFF(tempphaseim,strcat(g5folder,label_name,'_',num2str(regionidx),'_g5.tif'),'uint16');%Save the phase image
                       %  writeTIFF(labelmap,strcat(g5folder,label_name,'_lbl_',num2str(regionidx),'_g5.tif'));%Save the phase image
                       %writeTIFF(temptextonidx,strcat(g5folder,label_name,'_texidx_',num2str(regionidx),'_g5.tif'));%Save the phase image
                        save(strcat(g5folder,label_name,'_coord_',num2str(regionidx),'_g5.mat'),'r1','r2','c1','c2');
           
                   elseif (class==-1.0)
                       % writeTIFF(tempphaseim,strcat(nmfolder,label_name,'_small_',num2str(regionidx),'_nm.tif'),'uint16');%Save the phase image
                          
                       % writeTIFF(tempphaseim,strcat(nmfolder,label_name,'_',num2str(regionidx),'_nm.tif'),'uint16');%Save the phase image
                       %  writeTIFF(labelmap,strcat(nmfolder,label_name,'_lbl_',num2str(regionidx),'_nm.tif'));%Save the phase image
                       %  writeTIFF(temptextonidx,strcat(nmfolder,label_name,'_texidx_',num2str(regionidx),'_nm.tif'));%Save the phase image
                         save(strcat(nmfolder,label_name,'_coord_',num2str(regionidx),'_nm.mat'),'r1','r2','c1','c2');
           
                   elseif (class==-3.0) %BPH case
                       %  writeTIFF(tempphaseim,strcat(bphfolder,label_name,'_small_',num2str(regionidx),'_bph.tif'),'uint16');%Save the phase image
                       %  writeTIFF(tempphaseim,strcat(bphfolder,label_name,'_',num2str(regionidx),'_bph.tif'),'uint16');%Save the phase image
                       %  writeTIFF(labelmap,strcat(bphfolder,label_name,'_lbl_',num2str(regionidx),'_bph.tif'));%Save the phase image
                       %  writeTIFF(temptextonidx,strcat(bphfolder,label_name,'_texidx_',num2str(regionidx),'_bph.tif'));%Save the phase image
                       save(strcat(bphfolder,label_name,'_coord_',num2str(regionidx),'_bph.mat'),'r1','r2','c1','c2');
           
                   elseif (class==-5.0) %HGPIN
                        %writeTIFF(tempphaseim,strcat(hgpfolder,label_name,'_small_',num2str(regionidx),'_hgp.tif'),'uint16');%Save the phase image
                        %writeTIFF(tempphaseim,strcat(hgpfolder,label_name,'_',num2str(regionidx),'_hgp.tif'),'uint16');%Save the phase image
                        % writeTIFF(labelmap,strcat(hgpfolder,label_name,'_lbl_',num2str(regionidx),'_hgp.tif'));%Save the phase image
                        %writeTIFF(temptextonidx,strcat(hgpfolder,label_name,'_texidx_',num2str(regionidx),'_hgp.tif'));%Save the phase image
                         save(strcat(hgpfolder,label_name,'_coord_',num2str(regionidx),'_hgp.mat'),'r1','r2','c1','c2');
           
                   end
               end
          end
     else
         tempphaseim=0;         
     end
 
end