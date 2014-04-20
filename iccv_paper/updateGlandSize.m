function updateGlandSize(filenames,labelpath,glandsizepath)
    %This function computes the size of different glands in the biopsies
    %core...The glandsize as well as the distribution of glandsize will be
    %plotted...This code works on the 2048 x 2048 image. Otherwise, the
    %image has to be rescaled...
    %Compute the parameters of the biopsies using the groundtruth label...
    addpath(labelpath);
    addpath(glandsizepath);
    minsizeinpixel = 200;
    h = waitbar(0,'Calculating size of the glands progress...');
    nfeatures = 32;
    nsamplesaved = 0;
    feat = zeros(0,nfeatures);
    dataclass = zeros(0,1);
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
             cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            label_file_name = strcat(label_name,'_resized.mat');
            load(strcat(labelpath,label_file_name));
            glandsize_file_name = strcat(label_name,'_size.mat');
            disp(['Labeling: ' label_name ' ...']);
        
            if (~exist(glandsize_file_name,'file'))                
                %Create a gland mask...
                glandmap = lblim==1;
                glandmap = imfill(glandmap,'holes');
                %Remove too small regions
                glandmap = bwareaopen(glandmap,minsizeinpixel);
                
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
                %Plot the boundary of regions
                newim = zeros(size(glandmap));
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
                
                %Caculate the lumen parameters...
                lumenarea = zeros(nregions,1); 
                lumenmap = glandmap.*(lblim==0); %Get lumens within convex hull
                %Clear out small lumen area at the gland's boundary
                lumenmap = imopen(lumenmap,strel('disk',5));
                lumenmap = imfill(lumenmap,'holes');
                %Remove areas of lumen with too few pixels
                lumenmap = bwareaopen(lumenmap,200,8);
                
                %Compute the lumen area for each gland...
                pixidxperregion = regionprops(glandmap,'PixelIdxList');
                for glandregionidx = 1:nregions
                    lumenarea(glandregionidx,1) = sum(lumenmap(pixidxperregion(glandregionidx,1).PixelIdxList));
                    glandarea(glandregionidx,1)=areaarr(glandregionidx,1).Area;
                end
                lumentoglandratio = lumenarea./glandarea;
                
                %Save the calculated parameters.... for each biopsy
                %For each features, we save the mean value, median value,
                %mode and standard deviation value...
                %For gland area
                glandareamean = mean(glandarea);
                glandareamode = max(glandarea);
                glandareamed = median(glandarea);
                glandareastd = std(glandarea);
                %For gland perimeter
                glandpermean = mean(glandperi);
                glandpermode = max(glandperi);
                glandpermed = median(glandperi);
                glandperstd = std(glandperi);
                %For gland solidity
                glandsoldmean = mean(glandsold);
                glandsoldmode = max(glandsold);
                glandsoldmed = median(glandsold);
                glandsoldstd = std(glandsold);
                %For standard deviation of gland's radius...
                glandstddistmean = mean(glandstddistancearr);
                glandstddistmode = max(glandstddistancearr);
                glandstddistmed = median(glandstddistancearr);
                glandstddiststd = std(glandstddistancearr);
                %Max gland radius
                glandmaxdistmean = mean(glandmaxdistancearr);
                glandmaxdistmode = max(glandmaxdistancearr);
                glandmaxdistmed = median(glandmaxdistancearr);
                glandmaxdiststd = std(glandmaxdistancearr);
                %Gland compactness
                glandcompactmean = mean(glandcompactness);
                glandcompactmode = max(glandcompactness);
                glandcompactmed = median(glandcompactness);
                glandcompactstd = std(glandcompactness);
                %Lumen area
                lumenareamean = mean(lumenarea);
                lumenareamode = max(lumenarea);
                lumenareamed = median(lumenarea);
                lumenareastd = std(lumenarea);
                %Lumen to gland ratio
                lumenglandmean = mean(lumentoglandratio);
                lumenglandmode = max(lumentoglandratio);
                lumenglandmed = median(lumentoglandratio);
                lumenglandstd = std(lumentoglandratio);
                
                %Save the feature vector and the label...
                feat(end+1,:) = [glandareamean,glandareamode,glandareamed,glandareastd,...
                glandpermean,glandpermode,glandpermed,glandperstd,...
                glandsoldmean,glandsoldmode,glandsoldmed,glandsoldstd,...
                glandstddistmean,glandstddistmode,glandstddistmed,glandstddiststd,...
                glandmaxdistmean,glandmaxdistmode,glandmaxdistmed,glandmaxdiststd,...
                glandcompactmean,glandcompactmode,glandcompactmed,glandcompactstd,...
                lumenareamean,lumenareamode,lumenareamed,lumenareastd,...
                lumenglandmean,lumenglandmode,lumenglandmed,lumenglandstd...
                ];
                dataclass(end+1,:)= classidx;
                nsamplesaved =  nsamplesaved+1;
                figure(1);
                hold on;
                plot(log10(feat(end,:)));

                disp('Done.');
            end
       end
    end
    close(h);
    save('glandmorp.mat','feat','dataclass');
end