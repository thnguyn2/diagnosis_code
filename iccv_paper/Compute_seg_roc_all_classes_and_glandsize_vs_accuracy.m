%This source code perform multi-resolution segmentation for glands vs
%stroma. First, it will swipe through all folder, look for two fuzzy map of
%class prob at two different radii. Combine it and produce a final
%segmentation map
clc;
clear all;
close all;
datapath = '/Volumes/New_Athena/Dino_data/TMA_cores_and_diagnosis/';
texton_dir=strcat(datapath,'texdir/');
label_dir=strcat(datapath,'label/');

addpath(strcat(cd(cd('..')),'/support'));
[filenames,glandnames,classlist]=findFileNameFromROIs(datapath);
radius1 = 60;
radius2 = 90;
ngrades = size(filenames,1);

%The following parameters are good for diagnosis
glandthresh=0.43;
stromathresh=0.57;

%The following parameters are for the segmentation accuracy
%glandthresh=0.50;
%stromathresh=0.55;
fullcf = zeros(3,3);
colorarr = 'rbgykm';
gland_area_arr = zeros(0,1);
gland_accu_arr = zeros(0,1);
for classidx=1:ngrades %Go through different classes
               nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
               classname = classlist{classidx};
               curcf=zeros(3,3);
               allscore=zeros(0,1);
               alllabel = zeros(0,1);
               for sampleIdx=1:nsamples
                    cur_file_name = filenames{classidx,1}{sampleIdx,1};
                    dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                    slash_pos = strfind(cur_file_name,'/');
                    label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                    disp(['Working on ' label_name ', Classidx: ' num2str(classidx) ', Sample Idx: ' num2str(sampleIdx)])
                    seg_im_name_bf1 = strcat(cur_file_name(1:end-4),'_seg_ws',num2str(radius1),'_fuzzy.tif');  
                    seg_im_name_bf2 = strcat(cur_file_name(1:end-4),'_seg_ws',num2str(radius2),'_fuzzy.tif'); 
                    seg_im_name = strcat(cur_file_name(1:end-4),'_seg_multi_res_3072_strict.tif'); 
                    seg_im_name_gt = strcat(cur_file_name(1:end-4),'_seg_gt_3072.tif'); 
                    tx_hist_file_name = strcat(label_name,'_texton_hist_',num2str(radius2),'.mat'); %File path to load the label image
                    load(strcat(texton_dir,tx_hist_file_name),'lblim');
                   
                    if (exist(seg_im_name_bf1)&exist(seg_im_name_bf2)&(~exist(seg_im_name)))
                        fmap_1=imread(seg_im_name_bf1); %Read the first map
                        fmap_1 = bilateralFilter(fmap_1,'SigmaColor',0.1,'SigmaSpace',5);
                               
                        fmap_2=imread(seg_im_name_bf2); %Read the first map
                        fmap_2 = bilateralFilter(fmap_2,'SigmaColor',0.1,'SigmaSpace',40);
                     
                         lblim_bf=zeros(size(lblim)); %This is the new segmentation_map
                        nonlumenidx=find(lblim~=0);
                        glandidx=intersect(find(fmap_2<glandthresh),nonlumenidx); %Get the map for gland, make sure we are confident enough on the result

                        glandmap = zeros(size(lblim));
                        %Try to make sure that all glands we make are
                        %correct. It's better to under estimate it
                        %rather than overrestimate it.
                        glandmap(glandidx)=1;
                        glandmap =bwareaopen(glandmap,1000);
                        %Elminate very small hole
                        glandmapinv = 1-glandmap;
                        glandmapinv = bwareaopen(glandmapinv,200);
                        glandmap = 1-glandmapinv;
                        glandmap = imopen(glandmap,strel('disk',5));
                        glandidx = find(glandmap==1);
                        lblim_bf(glandidx)=1;
                        

                        stromaidx = intersect(find(fmap_2>stromathresh),nonlumenidx);
                        stromamap = zeros(size(lblim));
                        stromamap(stromaidx)=1;
                        stromamapinv=1-stromamap;
                        stromamapinv=bwareaopen(stromamapinv,800); %Get rid of small holes in stroma
                        stromamap = 1-stromamapinv;
                        stromaidx = find(stromamap==1);
                        lblim_bf(stromaidx)=2;
                        
                        non_det_idx = intersect(find(lblim~=0),find(lblim_bf==0));%Find all unsure region
                        lblim_bf(non_det_idx)=0;
                        %Determine more unsure pixels with the finer map,
                        %resolve small area of stroma
                        stromaidx2 = intersect(find(fmap_1>0.6),non_det_idx);
                        lblim_bf(stromaidx2)=2;
                        nonlumenidx = find(lblim_bf~=0);
                        curscore=fmap_2(nonlumenidx);
                        curlabelmap = lblim(nonlumenidx);
                        allscore(end+1:end+length(nonlumenidx))=curscore(:);
                        alllabel(end+1:end+length(nonlumenidx))=curlabelmap(:);
                        
                        figure(2)
                        subplot(121);imagesc(lblim_bf);title('Current multi-res segmentation');drawnow;
                        subplot(122);imagesc(lblim);title('Ground truth');drawnow;
                        writeTIFF(lblim_bf,seg_im_name);
                        writeTIFF(lblim,seg_im_name_gt);
                    
                        
                    
                    end
                    
                    if (exist(seg_im_name_gt,'file')&exist(seg_im_name,'file'))
                    
                        lblim = imread(seg_im_name_gt);
                        lblim_bf = imread(seg_im_name);

                        [cm,gorder] = confusionmat(lblim(:),lblim_bf(:));
                        curcf = curcf+cm;  
                        curcf_disp = curcf./repmat(sum(curcf,2),1,3);                   
                        disp(['Total number of pixels: ' num2str(sum(sum(curcf)))])
                        figure(1);
                        imagesc(curcf_disp);title('Current confusion matrix'); colorbar;drawnow;

                        %Compute grandsize vs accuracy
                        figure(1);
                        subplot(121);imagesc(lblim);
                        subplot(122);imagesc(lblim_bf);
                        glandmap = imfill(lblim==1,'hole');
                        glandpixidxlist = regionprops(glandmap,'PixelIdxList');
                        nglands_org = size(glandpixidxlist,1);
                        for glandidx = 1:nglands_org
                           %Get the number of pixels in each glands
                           pixlist=glandpixidxlist(glandidx).PixelIdxList;
                           npix = length(pixlist);
                           if (npix>20000)
                                cur_accu = sum(lblim(pixlist)==lblim_bf(pixlist))/npix;
                                gland_area_arr(end+1)=npix;
                                gland_accu_arr(end+1)=cur_accu;
                           end
                        end
                        min_area = min(gland_area_arr);
                        max_area = 0.15*max(gland_area_arr);
                        nbins = 50;
                        area_arr = linspace(min_area,max_area,nbins);
                        [area_hist,area_coord]=hist(gland_area_arr,area_arr);
                        area_hist = area_hist/max(area_hist(:));
                        gland_accu_hist = zeros(1,nbins);
                        %Compute the bin_idx for each datasample
                        bin_area_size = (max_area-min_area)/nbins;
                        gland_bin_idx = min(round((gland_area_arr-min_area)/bin_area_size+0.5),nbins);
                        for binidx=1:nbins
                            curidx = find(gland_bin_idx==binidx);
                            if (~isempty(curidx))
                                gland_accu_hist(binidx)=mean(gland_accu_arr(curidx));
                            end
                        end
%                         [sorted_area,I] = sort(gland_area_arr,'descend');
%                            sorted_accu = gland_accu_arr(I);
%                            %Draw the frequency of the area
%                            
%                            [hist_val,area_coord]=hist(sorted_area,nbins);
%                            hist_val = hist_val/max(hist_val);
                            figure(2);
                            plot(area_coord,area_hist,'-b','linewidth',2);hold on;
%                            nglandsperinterval = ceil(length(sorted_accu)/nbins);
%                            mean_accu_hist =imresize(sorted_accu,[1 nbins]);
                            figure(2);
                            plot(area_coord,gland_accu_hist,'-r','linewidth',2);hold off;
                            legend('Area frequency','Accuracy');
                            drawnow;
                    end
               end
               curcf_disp
%                nsamplesmax = min(100000,length(alllabel(:)));
%                p = randperm(length(alllabel(:)));
%                alllabel = alllabel(p(1:nsamplesmax));
%                allscore = allscore(p(1:nsamplesmax));
%                
%                [xs,ys,t,auc]=perfcurve(alllabel,allscore,2);
%                figure(3);
%                plot(xs,ys,colorarr(mod(classidx,6)+1));
%                xlabel('Threshold');
%                xlabel('P_F (stroma)'); ylabel('P_D (gland)');
                

               %save(strcat('./confusionmatrices/',classname,'cf.mat'),'curcf','curcf_disp','xs','ys','t','auc');
%               fullcf = fullcf + curcf;

               %
               
end
save('area_vs_accuracy.mat','area_coord','area_hist','gland_accu_hist','gland_area_arr','gland_accu_arr');
%  fullcf_disp = fullcf./repmat(sum(fullcf,2),1,3);
  %save(strcat('.\confusionmatrices\fullcf.mat'),'fullcf','fullcf_disp');
           
                 
