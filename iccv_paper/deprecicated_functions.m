function [feat,dataclass,samplename]=mean_size_analysis(filenames,labelpath,mode,gradelist)
      feat = zeros(0,7);
      dataclass = cell(0,1); %This is the class of the sample.
      ngrades = size(filenames,1);
      samplename = cell(0,1);
      minglandsizeinpixel = 2000; %Define the minimum size of the gland in pixels
      totalnumberofsampleprocessed = 0;
      for classidx=ngrades:-1:1 %Go through different classes
           nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
           for sampleIdx=1:min(nsamples,nsamples)
                cur_file_name = filenames{classidx,1}{sampleIdx,1};
                %Check to see if the label file exist
                dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                slash_pos = strfind(cur_file_name,'\');
                label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                samplename{end+1,:}=label_name;
                label_file_name = strcat(label_name,'_resized.mat');
                disp(['Working on : ' label_name ' ...']);
                smallim_name = strcat(cur_file_name(1:end-4),'_small.tif');
                if (strcmp(mode,'gt')==1)
                    load(strcat(labelpath,label_file_name));
                else
                    lblim = imread(strcat(cur_file_name(1:end-4),'_seg_multi_res_3072.tif'));                
                end
                figure(3);
                imagesc(lblim);
                drawnow;
                %Create a gland mask...
                glandmap = lblim==1;
                glandmap = imfill(glandmap,'holes');
                %Remove too small regions
                glandmap = bwareaopen(glandmap,minglandsizeinpixel);

                %Caculate parameters for glands...
                areaarr = regionprops(glandmap,'Area'); %1
                perimeterarr = regionprops(glandmap,'Perimeter'); %3

                solidityarr = regionprops(glandmap,'Solidity');%2
                %Next, compute the distance ratio, which is defined as the
                %ratio of the distance between the max distance from the
                %centroid
                %First, get the controid coordinates
                centroidarr = regionprops(glandmap,'Centroid');
                %Get boundary pixels
                boundaryarr = bwboundaries(glandmap,'noholes');
                nregions = size(centroidarr,1); %Number of centroid regions

                glandmaxdistancearr = zeros(nregions,1); %4
                glandmindistancearr = zeros(nregions,1); 
                glandstddistancearr = zeros(nregions,1); %5
                glandsold           = zeros(nregions,1);
                glandperi           = zeros(nregions,1);
                glandarea           =zeros(nregions,1); 

                %Compute the compactness of different glands...
                glandcompactness = zeros(nregions,1);
                for regionidx=1:nregions
                    cur_cent = (centroidarr(regionidx,1).Centroid);    %First column is X, second column is Y          
                    curbound2d =  boundaryarr{regionidx,1}; %First column is Y, second column is X
                    distance = repmat([cur_cent(2) cur_cent(1)],[size(curbound2d,1) 1])-curbound2d;
                    distance = sqrt(sum(distance.^2,2));
                    glandmaxdistancearr(regionidx,1) = max(distance)/mean(distance);
                    glandmindistancearr(regionidx,1) = min(distance)/mean(distance);
                    glandstddistancearr(regionidx,1) = std(distance)/mean(distance);
                    glandcompactness(regionidx,1) = (size(curbound2d,1)^2)/areaarr(regionidx,1).Area;%6
                    glandsold(regionidx,1)=solidityarr(regionidx,1).Solidity;
                    glandperi(regionidx,1)=perimeterarr(regionidx,1).Perimeter;
                    glandarea(regionidx,1) = areaarr(regionidx,1).Area;
                end

                %Calcuate the lumen parameters...
                lumenarea = zeros(nregions,1); 
                nonsamplemap = (lblim~=0);
                nonsamplemap = imclose(nonsamplemap,strel('disk',30)); %Avoid the phenomenon when the lumen is cut through
                nonsamplemap = imfill(nonsamplemap,'holes');
                
                
                lumenmap = nonsamplemap.*(lblim==0); %Get lumens within convex hull
                %Clear out small lumen area at the gland's boundary
                lumenmap = imopen(lumenmap,strel('disk',10));
                lumenmap = imfill(lumenmap,'holes');
                %Remove areas of lumen with too few pixels
                lumenmap = bwareaopen(lumenmap,2000,8);

                %Compute the lumen area for each gland...
                pixidxperregion = regionprops(glandmap,'PixelIdxList');
                for glandregionidx = 1:nregions
                     lumenarea(glandregionidx,1) = sum(lumenmap(pixidxperregion(glandregionidx,1).PixelIdxList));
                     glandarea(glandregionidx,1)=areaarr(glandregionidx,1).Area;
                end
                lumentoglandratioarr = lumenarea./glandarea;
                %For gland area
                glandareamean = mean(glandarea);
                %For gland perimeter
                glandpermean = mean(glandperi);
                %For gland solidity
                glandsoldmean = mean(glandsold.*glandarea)/glandareamean; %Use a weighted mean instead of the normal mean
                %For standard deviation of gland's radius...
                glandstddistmean = mean(glandstddistancearr.*glandarea)/glandareamean;
                %Max gland radius
                glandmaxdistmean = mean(glandmaxdistancearr.*glandarea)/glandareamean;
                %Gland compactness
                glandcompactmean = mean(glandcompactness.*glandarea)/glandareamean;
                lumentoglandmean = mean(lumentoglandratioarr.*glandarea)/glandareamean;
                
      
                disp(['Class - #glands - Area - Perimeter - Solidity - Distance std - Max distance - Compactness - Lumen/area']);
                disp([gradelist(classidx) num2str(nregions) num2str(glandareamean) num2str(glandpermean)  num2str(glandsoldmean)...
                      num2str(glandstddistmean) num2str(glandmaxdistmean) num2str(glandcompactmean) num2str(lumentoglandmean)]);

                %Save the feature vector and the label...
                feat(end+1,:) = [glandareamean,glandpermean,glandsoldmean,...
                    glandstddistmean,glandmaxdistmean,glandcompactmean,...
                    lumentoglandmean];
                dataclass(end+1,:)= gradelist(classidx);
                
                %Plot the feature distribution with a colorbar
                [meanNM,meanBPH,meanHGP,meanHGC,meanLGC,stdNM,stdBPH,stdHGP,stdHGC,stdLGC]=compute_feature_statistics(feat,dataclass);
                %Normalize w.r.t all classes so that the maximum value is 1
                mean_arr = [meanNM;meanBPH;;meanLGC;meanHGP;meanHGC];
                std_arr = [stdNM;stdBPH;stdLGC;stdHGP;stdHGC];
                figure(2);
                barwitherr(std_arr', mean_arr');
                legend('NM','BPH','LGC','HGPIN','HGC');
                set(gca,'XTickLabel',{'Area','Peri.','Sold.','DistStd','DistMax','Compact.','Lumen/Gland'})
                drawnow;
                totalnumberofsampleprocessed = totalnumberofsampleprocessed + 1;
                disp(['Done...' num2str(totalnumberofsampleprocessed) ' samples']);
            
              
           end
      end      
    
end
