function truncate_segmented_map_with_from_core_label
    clc;
    clear all;
    close all;
    datapath = 'E:\TMA_cores_and_diagnosis\';
    texton_dir=strcat(datapath,'texdir\');
    label_dir=strcat(datapath,'label\3class_results\');
    core_groundtruth_dir = strcat(datapath,'diagnosis_of_vicky\');%This is the groundtruth for the core
    addpath('.\support')
    addpath(strcat(cd(cd('..')),'\support'));
    [filenames,glandnames,classlist,alltiffilelist]=findFileNameFromROIs(datapath);
    noglandmarkuplist = {'N1','O1','P1','S1','J2','N2','P2','S2','I3','K3','M3','N3','P3','T3',...
        'Y3','A4','I4','Q4','S4','T4','Y4','A6','I6','Q6','G7','I7','K7','M7','N7','O7','T7',...
        'H8','M8','N8','U8','W8','A9','I9','N9','T9','W9','D11','I11','M11','N11','P11','T11',...
        'Z11','D12','P12','Q12','U12','W12','Z12','D13','K13','M13','Q13','T13','W13','Z13',....
        'D14','J14','K14','M14','Q14','T14','W14','Z14','Z15','K16','T16','W16','Y16','Z16',...
        'C17','S17','Y17','Z17','U18','W18','G19','I19','U19','W19','Y19','Z19','AA19'};%This is the list of the cores that doesn't have the student marking
    notusenonegtdata = 1;
    for fileidx = 1:length(alltiffilelist)
        %cur_file_name = alltiffilelist{fileidx};
        cur_file_name = 'E:\TMA_cores_and_diagnosis\dHGPIN\AB11.tif';
        dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
        slash_pos = strfind(cur_file_name,'\');
        label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
        if ((notusenonegtdata)&&(~ismember(label_name,noglandmarkuplist)))

            disp(['Working on ' label_name ', Sample Idx: ' num2str(fileidx)])
            phaseim = imread(strcat(cur_file_name(1:end-4),'_small.tif'));%Downsampled image  
            labelim3cl = imread(strcat(label_dir,label_name,'_90_3cl.tif'));%3 class map - The segmentation map needs to be clear
            %%more on cleaning the segmentation results. The lumen is perfect
            %%(March 18, 15)
            labelim2cl = imread(strcat(cur_file_name(1:end-4),'_seg_multi_res_3072_strict_hf.tif')); %2 class map.

            %Make sure that there is no holes inside the glands and inside the
            %stroma. This code should be revised so that only a 3 class
            %diagnosis results can be used
            lumenmap = (labelim3cl==0);
            stromamap = (labelim2cl==2.0);
            %Clear out small regions
            invstromamap = 1-stromamap;
            invstromamap = bwareaopen(invstromamap,5000);
            stromamap = 1-invstromamap;
            glandmap = (labelim2cl==1.0);
            glandmap = imfill(glandmap,'holes');
            glandidx = find(glandmap == 1);
            labelim = 2*stromamap;
            labelim(glandidx)=1;
            labelim = labelim .*(1-lumenmap); %Mask out lumen puxe
            %Load the core diagnosis and use it to crop the segmented map
            coregt_map_file_name = strcat(core_groundtruth_dir,label_name,'_roi_diag.tif');
            core_index_texton_name = strcat(texton_dir,label_name,'_texton_index_map.tif');
            if (exist(coregt_map_file_name))
                coregt_map = imread(coregt_map_file_name);
                figure(1);
                subplot(121);
                imagesc(coregt_map);colorbar
                subplot(122);
                imagesc(labelim);colorbar
                texidxmap = imread(core_index_texton_name);%Load the histogram of texton



                drawnow
                create_crop_imagesc(3,coregt_map,phaseim,texidxmap,labelim,label_name,datapath);%Create all the image of grade 3 cores
                create_crop_imagesc(4,coregt_map,phaseim,texidxmap,labelim,label_name,datapath);%Grade 4 case
                create_crop_imagesc(5,coregt_map,phaseim,texidxmap,labelim,label_name,datapath);%Grade 5 case
                create_crop_imagesc(-1,coregt_map,phaseim,texidxmap,labelim,label_name,datapath);%Normal case
                create_crop_imagesc(-3,coregt_map,phaseim,texidxmap,labelim,label_name,datapath);%BPH case


            else
                disp(['Core ' label_name ' ground truth does not exist...'])
            end
        else
            disp(['Skip ' label_name 'for belonging to the exclusion list'])
        end

    end
end

function [tempphaseim,temptextonidx]=create_crop_imagesc(class,coregt_map,phaseim,texidxmap,labelim,label_name,datapath)
     g3folder = strcat(datapath,'diagnosis_of_vicky\g3\');
     g4folder = strcat(datapath,'diagnosis_of_vicky\g4\');
     g5folder = strcat(datapath,'diagnosis_of_vicky\g5\');
     nmfolder = strcat(datapath,'diagnosis_of_vicky\nm\');
     bphfolder = strcat(datapath,'diagnosis_of_vicky\bph\');
    
     bw = (coregt_map==class);
     if (sum(bw(:))>0)
          regions = regionprops(bw,'PixelIdxList'); %Get all the grade 3 regions
          areas = regionprops(bw,'Area');
          for regionidx=1:length(regions)
               curidx1 = regions(regionidx).PixelIdxList; 
               curarea = areas(regionidx).Area;
               if (curarea>=30000)
                   [currowidx,curcolidx]=ind2sub(size(coregt_map),curidx1);
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
                   tempphaseim = zeros(bbnrows,bbncols);
                   tempphaseim(curidx)=phaseim(curidx1);
                   labelmap(curidx) = labelim(curidx1);
                   temptextonidx(curidx) = texidxmap(curidx1);
                   figure(3);
                   imagesc(tempphaseim);drawnow;
                   figure(4);
                   imagesc(labelmap);drawnow;
                   figure(5);
                   imagesc(temptextonidx);drawnow;
                   if (class==4)
                         writeTIFF(tempphaseim,strcat(g4folder,label_name,'_small_',num2str(regionidx),'_g4.tif'));%Save the phase image
                         writeTIFF(labelmap,strcat(g4folder,label_name,'_lbl_',num2str(regionidx),'_g4.tif'));%Save the phase image
                         writeTIFF(temptextonidx,strcat(g4folder,label_name,'_texidx_',num2str(regionidx),'_g4.tif'));%Save the phase image
                  
                   elseif (class==3)
                         writeTIFF(tempphaseim,strcat(g3folder,label_name,'_small_',num2str(regionidx),'_g3.tif'));%Save the phase image
                         writeTIFF(labelmap,strcat(g3folder,label_name,'_lbl_',num2str(regionidx),'_g3.tif'));%Save the phase image
                         writeTIFF(temptextonidx,strcat(g3folder,label_name,'_texidx_',num2str(regionidx),'_g3.tif'));%Save the phase image
                  
                   elseif (class==5)
                         writeTIFF(tempphaseim,strcat(g5folder,label_name,'_small_',num2str(regionidx),'_g5.tif'));%Save the phase image
                         writeTIFF(labelmap,strcat(g5folder,label_name,'_lbl_',num2str(regionidx),'_g5.tif'));%Save the phase image
                         writeTIFF(temptextonidx,strcat(g5folder,label_name,'_texidx_',num2str(regionidx),'_g5.tif'));%Save the phase image
                  
                   elseif (class==-1.0)
                         writeTIFF(tempphaseim,strcat(nmfolder,label_name,'_small_',num2str(regionidx),'_nm.tif'));%Save the phase image
                         writeTIFF(labelmap,strcat(nmfolder,label_name,'_lbl_',num2str(regionidx),'_nm.tif'));%Save the phase image
                         writeTIFF(temptextonidx,strcat(nmfolder,label_name,'_texidx_',num2str(regionidx),'_nm.tif'));%Save the phase image
                  
                   elseif (class==-3.0) %BPH case
                         writeTIFF(tempphaseim,strcat(bphfolder,label_name,'_small_',num2str(regionidx),'_bph.tif'));%Save the phase image
                         writeTIFF(labelmap,strcat(bphfolder,label_name,'_lbl_',num2str(regionidx),'_bph.tif'));%Save the phase image
                         writeTIFF(temptextonidx,strcat(bphfolder,label_name,'_texidx_',num2str(regionidx),'_bph.tif'));%Save the phase image
                  
                   end
               end
          end
     else
         tempphaseim=0;         
      end
end