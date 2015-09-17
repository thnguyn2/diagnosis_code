function generate_gland_image_from_roi
    %This function reads the ROIs and the phase image. Then, it generates
    %different gland images. Note that this function only support
    %Freehand ROIS
    
    %First, look for all the files of different classes from the list of
    %available rois
    clc;
    clear all;
    close all;
    allcorepath='H:\TMA_cores_and_diagnosis\';
    roipath = strcat(allcorepath,'Verified_rois\Verified_diagnosis_results\');
    %First, look for the list of all cores
     %Look for all file names
    filenames =  dir(strcat(roipath,'*.zip'));
    nfiles = length(filenames);
    filelist = cell(0,1);
    for filenameidx=1:nfiles
        filelist{end+1,1} =strcat(roipath,filenames(filenameidx,1).name);
    end
    
    %Step 2: Search for all tif files in all folders
    tiffilelist = cell(0,1);
    tiffilelistnopath = cell(0,1);
    folderlist = dir(strcat(allcorepath,'d*'));%List all subfolders in the current folder
    nfolders = size(folderlist,1);%Number of subfolder
    for folderidx=1:nfolders
        curfoldername = strcat(allcorepath,folderlist(folderidx,1).name,'\');
        %Search all the tif files in the current folder
        tiffilenames = dir(strcat(curfoldername,'*.tif'));
        for tiffileidx=1:size(tiffilenames,1)
            tiffilelistnopath{end+1,1}=tiffilenames(tiffileidx,1).name;
            tiffilelist{end+1,1} =strcat(curfoldername,tiffilenames(tiffileidx,1).name);
        end        
    end
    
    ntiffiles = size(tiffilelist,1);
    h = waitbar(0,'Processing data...');
    for fileIdx = 1:nfiles
        waitbar(fileIdx/nfiles,h,'Progress...')
        cur_file_name = filelist{fileIdx,1};
        waitbar(fileIdx/nfiles,h,'Progress...')
        cur_file_name = 'H:\TMA_cores_and_diagnosis\Verified_rois\Verified_diagnosis_results\33A8.zip';
        dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
        slash_pos = strfind(cur_file_name,'\');
        label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
        disp(['Working on ' label_name]);
        label_name = label_name(3:end); %Truncate the first 2 characters, which are for Gleason grading
        %Generate the tif file name to search
        tif_name = strcat(label_name,'.tif');
        %Search for match tif file
        foundidx=-1;
        for tiffileidx=1:ntiffiles
            if (strcmp(tif_name,tiffilelistnopath{tiffileidx,1}))
                foundidx=tiffileidx;
                break;
            end
        end
     
        
        %Step 3: Red the roi and the tiff file
        roi = ReadImageJROI(filelist(fileIdx,1));
        nroi = size(roi{1,1},2);
        if (foundidx~=-1)
            curtifpath = tiffilelist{foundidx,1};
            glandimfullpath = strcat(curtifpath(1:end-4),'Glands.tif');
            errorcount = 0;
            if (~exist(glandimfullpath))
%                 curtifim = imread(curtifpath);
%                 nrows = size(curtifim,1);
%                 ncols = size(curtifim,2);
% 
%                 %Generate the mask for ROI
%                 mask = zeros(nrows,ncols);
%                 glandim =ones(nrows,ncols)*8192;
                for roiidx=1:nroi
                    if (strcmp(roi{1,1}{1,roiidx}.strType,'Freehand'))
                        curbound = roi{1,1}{1,roiidx}.mnCoordinates; %The columns are [X, Y]
%                         BW = roipoly(mask,curbound(:,1),curbound(:,2));
%                         mask = mask + BW;
                    else
                        errorcount = errorcount + 1;
                        disp(['[Error]: ',roi{1,1}{1,roiidx}.strName]);
                    end

                end
%                 glandpix = find(mask~=0);
%                 glandim(glandpix) = curtifim(glandpix);
% 
%                 glandim = cast(glandim,'uint16');
%                 %Save the gland image in the same folder
%                 disp('...Saving...');
% 
%                 figure(1);
%                 subplot(121);
%                 imagesc(curtifim);
%                 subplot(122);
%                 imagesc(glandim);
% 
%                 writeTIFF(glandim, glandimfullpath,'uint16')
            else
                for roiidx=1:nroi
                    if (~strcmp(roi{1,1}{1,roiidx}.strType,'Freehand'))
                        errorcount = errorcount + 1;
                    end                
                end
            end
        else
            disp(strcat('*********Cannot find: ',tif_name,'******'));
        end
        disp(['Error count: ' num2str(errorcount) '/' num2str(nroi)]);
    end
    close(h);
    
end