%This source code performs the region growing segmentation in order to get
%all glands precisely....
%segmentation map
clc;
clear all;
close all;
datapath = 'H:\TMA_cores_and_diagnosis\';
texton_dir=strcat(datapath,'texdir\');
label_dir=strcat(datapath,'label\2class_results\');

addpath(strcat(cd(cd('..')),'\support'));
[filenames,glandnames,classlist,alltiffilelist]=findFileNameFromROIs(datapath);
radius1 = 60;
radius2 = 90;
ngrades = size(filenames,1);

%The following parameters are good for diagnosis
glandthresh=0.3;%Make sure that all glands we get are correct
stromathresh=0.5;

%The following parameters are for the segmentation accuracy
%glandthresh=0.50;
%stromathresh=0.55;
fullcf = zeros(3,3);
colorarr = 'rbgykm';
retrain = 1;
fillinlumeninsideglandandstroma = 1
for fileidx = 1:length(alltiffilelist)
                    cur_file_name = alltiffilelist{fileidx};
                    dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                    slash_pos = strfind(cur_file_name,'\');
                    %label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                    label_name = 'O2';
                    disp(['Working on ' label_name ', Sample Idx: ' num2str(fileidx)])
                    seg_im_name_bf1 = strcat(label_dir,label_name,'_',num2str(radius2),'_2cl.tif');  
                    seg_im_name = strcat(label_dir,label_name,'_seg_rg.tif'); %Save name for the region growed image

                    if ((exist(seg_im_name_bf1)&(~exist(seg_im_name)))|(retrain==1))
                        fmap_1=imread(seg_im_name_bf1); %Read the first map
%                         h = fspecial('gaussian',50,20)
%                         fmap_2 = imfilter(1-fmap_1,h,'same');
%                         lumenmap = fmap_1==0;
%                         se = strel('disk', 30);
%                         fmap_4 = imclose(fmap_2, se);
%                         fmap_4 = imregionalmax(fmap_4);
%                         se2 = strel(ones(5,5));
%                         fmap_4 = imopen(fmap_4, se2);
%                         coreregion = imfill(1-lumenmap,'holes');%Get rid of all holes
%                         coreregion = imfill(imclose(coreregion,strel('disk',30)),'holes');
%                         fmap_4 =fmap_4.*coreregion;
%                         fmap_4 = bwareaopen(fmap_4,5000);
%                         figure
%                         imagesc(fmap_4);
%                         figure
%                         imshow(Io), title('Opening (Io)')
%                         h = fspecial('gaussian',100,50)
%                         fmap_3 = imfilter(fmap_1,h,'same');
%                         lumenmap = (fmap_1==0); %find all the lumenpixel
%                         lumenmap = bwareaopen(lumenmap,500);
%                         fmap_2 = fmap_1;
%                        fmap_2(find(lumenmap==1))=3;
                        lblim_bf=zeros(size(fmap_1)); %This is the new segmentation_map                                           
                        glandmap=fmap_1<glandthresh;
                        glandmap =bwareaopen(glandmap,1000);
                        %Elminate very small hole inside the gland
                        glandmapinv = 1-glandmap;
                        glandmapinv = bwareaopen(glandmapinv,3000);%Eliminate all the holes in glands. Set this number larger will eliminate larger hole
                        glandmap = 1-glandmapinv;
                        glandmap = imopen(glandmap,strel('disk',15));%Make cut that we cut out all the connection between overlapping
                        glandidx = find(glandmap==1);
                        lblim_bf(glandidx)=1;
                        %Remove very small gland area
                        lblim_open = bwareaopen(lblim_bf,5000);
                        lblim_bf = lblim_open + 2*(lblim_bf-lblim_open);%Set the difference pixels to stroma
                        
                        %Now, work on stroma region
                        stromaidx = intersect(find(fmap_2>stromathresh-0.1),nonlumenidx);
                        stromamap = zeros(size(fmap_1));
                        stromamap(stromaidx)=1;
                        stromamap = imdilate(stromamap,strel('disk',10));
                        stromamapinv=1-stromamap;
                        stromamapinv=bwareaopen(stromamapinv,2000); %Get rid of small holes in stroma
                        stromamap = 1-stromamapinv;
                        stromaidx = find(stromamap==1);
                        lblim_bf(stromaidx)=2;
                        
%                         non_det_idx = intersect(find(fmap_1~=0),find(lblim_bf==0));%Find all unsure region
%                         stromaidx2 = intersect(find(fmap_1>0.5),non_det_idx);
%                         lblim_bf(stromaidx2)=2;
                        
                        figure(1);
                        subplot(121);imagesc(fmap_1);title('Scale 1');
                        subplot(122);imagesc(fmap_2);title('Scale 2');drawnow
                        figure(2)
                        imagesc(lblim_bf);title('Current multi-res segmentation');drawnow;
                        nonlumenidx = find(lblim_bf~=0);
                        writeTIFF(lblim_bf,seg_im_name);

                    end
               
   end
%  fullcf_disp = fullcf./repmat(sum(fullcf,2),1,3);
  %save(strcat('.\confusionmatrices\fullcf.mat'),'fullcf','fullcf_disp');
           
                 
