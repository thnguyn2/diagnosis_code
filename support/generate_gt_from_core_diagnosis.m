clc;
clear all;
close all;
disp('Generating groundth truth data from core diagnosis results...')
disp('Reading Excel label data...')
label_file_name = 'C:\Users\Tan\Dropbox\H & E for Vicky\label.xlsx';
[num_data,text_data,raw_data] = xlsread(label_file_name);
disp('Done.');
labeldir = 'G:\TMA_cores_and_diagnosis\diagnosis_of_vicky\';

%First, extract the core's name and the regions within the cores
ncores=length(raw_data(2:end,1));
disp(['Found ' num2str(ncores) ' different cores']);
corelist = raw_data(2:end,1);
coreregions = cell(ncores,1);
for coreidx=1:ncores
    rawcurlabel = raw_data(1+coreidx,2:end);
    tempdata=cell(0,1);
    for regionidx=1:length(rawcurlabel)%Find all the regions that is not NaN
        if (~isnan(rawcurlabel{regionidx}))
            tempdata{end+1,1}=rawcurlabel{regionidx};
        end
    end
    coreregions{coreidx,1}=tempdata;
end

%Step 2: extract the list of cores in the folder and create the cores' mask
corefolder = 'G:\TMA2D_all_cores\';
roifolder = 'G:\TMA_cores_and_diagnosis\diagnosis_of_vicky\'
tifffiles = dir(strcat(corefolder,'*.tif')); %Get all the tif files
tifflist = cell(0,1);
roilist = cell(0,1);
disp('Searching for all *.zip files...')
for tiffidx = 1:length(tifffiles(:,1))
    tifflist{end+1,1}=tifffiles(tiffidx).name(1:end-4);
    roilist{end+1,1}=strcat(roifolder,tifflist{end,1},'.zip');        
end
disp('Done.');

ncores = length(tifflist);
%[temp1,temp2,matchingidx]=intersect(tifflist,corelist);%Matchingidx is the array of indices in the corelist

%Step 3: read the tifffiles, downsample it and prepare the label map
nrows = 3072;
ncols = 3072;
curnrows=10000;
curncols=10000;
addpath('C:\Users\Tan\Dropbox\current_working_code\cancer_diagnosis\support\');%Add the support path to read imagej's roi
for tiffidx=1:ncores
   %Get the roi name
   curfilename = tifflist{tiffidx};
   disp(['Working on ' curfilename '...'])
   roiname = strcat(roifolder,curfilename,'.zip');%Check if we have the tif file and the roi file or not
   tifffilename = strcat(corefolder,curfilename,'.tif');
   if ((exist(roiname))&&(exist(tifffilename))&&(~exist(strcat(labeldir,curfilename,'_roi_diag.tif'))))
         sROI = ReadImageJROI(roiname);%Read the roi file
         %curim = imread(tifffilename);%Read the tiff file
         %[curnrows,curncols] = size(curim);
         %phaseresizedim=imresize(curim,[nrows ncols]);%Resize the image to match the new size
         nrois = length(sROI); %Number of regions marked
         newmap = zeros(nrows,ncols);
         roiareaslist = zeros(0,1); %This is the areas of the ROIs
         roiindexlist = cell(0,1);
         roitypelist = zeros(0,1);
         for roiIdx=1:nrois
            if (strcmp(sROI{1,roiIdx}.strType,'Freehand')) %just care about free hand selection
                xycoord = sROI{1,roiIdx}.mnCoordinates;
                xycoord(:,1)=round(xycoord(:,1)*ncols/curncols); %Coordinates of columns
                %Check for the boundary
                tempidx = find(xycoord(:,1)==0);
                xycoord(tempidx,1)=1;
                tempidx = find(xycoord(:,1)>=ncols);
                xycoord(tempidx,1)=ncols;
               
                xycoord(:,2)=round(xycoord(:,2)*nrows/curnrows);
                tempidx = find(xycoord(:,2)==0);
                xycoord(tempidx,2)=1;
                tempidx = find(xycoord(:,2)>=nrows);
                xycoord(tempidx,2)=nrows;
                
                
                bwmask = poly2mask(xycoord(:,1),xycoord(:,2),nrows,ncols);%Create a mask with specified boundary
                matchidx = find(strcmp(corelist,curfilename)==1);%Go look up to see what is the index of the matched core labeled
                %Further processing to find the type of rois here
                if ((nrois~=size(coreregions{matchidx,1},1)))
                    disp(['Number of cores marked: ' num2str(nrois)]);
                    disp(['Number of cores labeled in the excel files: ' num2str(size(coreregions{matchidx,1},1))]);
                    error('Number of cores marks and number of label does not match. Please double check')
                end
                inneridx = find(bwmask==1);
                roiareaslist(end+1) = length(inneridx);
                roiindexlist{end+1}=inneridx; %Add the current roi to the indexing list
                %Now, save the type of the regions as well
                curtype=coreregions{matchidx,1}{roiIdx,1};
                switch curtype
                    case {3}
                        roitype = 3.0;
                    case {34}
                        roitype = 3.4;
                    case {35}
                        roitype = 3.5;
                    case {4}
                        roitype = 4.0;
                    case {43}
                        roitype = 4.3;
                    case {45}
                        roitype = 4.5;
                    case {5}
                        roitype = 5.0;
                    case {54}
                        roitype = 5.4;
                    case {'NECRO'}
                        roitype = 6.0;
                    case {'SCRT'}
                        roitype = 7.0;
                    case {'CORPORA'}
                        roitype = 8.0;
                    case {'ATPA'}
                        roitype = 9.0;
                    case {'BVSL'}
                        roitype = 10.0;
                    case {'INFL'}
                        roitype = 11.0;
                    case {'NERV'}
                        roitype = 12.0;
                    case {'SQUAMOUS'}
                        roitype = 13;
                    case {'TUMOR'}
                        roitype = 3.0;
                    case {'NM'}
                        roitype = -1.0;
                    case {'ATPH'}
                        roitype = -2.0;
                    case {'BPH'}
                        roitype = -3.0;
                    case {'LGP'}
                        roitype = -4.0;
                    case {'HGP'}
                        roitype = -5.0;
                    case {'STROMA'}
                        roitype = -6.0;
                    case {'?'}
                        roitype = -10.0;
                   
                    otherwise
                        curtype
                        error('Unspecified type');
                        
                end
                roitypelist(end+1)=roitype;
            end
         end
         
         %Now, start arranging the size of the roi
         [temp,sortedidxlist]=sort(roiareaslist(:),'descend');
         %Now, start placing the roi in
         for sortedidx = 1:length(roiareaslist)
            oldidx = sortedidxlist(sortedidx);
            newmap(roiindexlist{oldidx})=roitypelist(oldidx);
         end
         
         figure(1);
         subplot(121);imagesc(newmap);colorbar
         %subplot(122);imagesc(phaseresizedim);
         drawnow
         writeTIFF(newmap, strcat(labeldir,curfilename,'_roi_diag.tif'));
       disp('Done.')
   else
       disp('Did not find the roi file...')
   end
   
   
end
