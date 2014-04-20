function updateGlandSizeNoGT(filenames,fimpath,morppath)
    %This function computes the size of different glands in the biopsies
    %core...The glandsize as well as the distribution of glandsize will be
    %plotted...This code works on the 2048 x 2048 image. Otherwise, the
    %image has to be rescaled...
    %The procedure done in the code will be as follows
    %Read the raw fim image => Calculate a label image => Calculate the
    %parameters...
    %Compute the parameters of the biopsies using the groundtruth label...
    addpath(fimpath);
    addpath(morppath);
    minsizeinpixel = 200;
    h = waitbar(0,'Calculating size of the glands progress...');
    nfeatures = 32;
    nsamplesaved = 0;
    feat = zeros(1,nfeatures);
    dataclass = zeros(0,1);
    totalfeat = zeros(0,nfeatures);
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            label_file_name = strcat(label_name,'_svm_200k_ws40_1024.mat');
            load(strcat(fimpath,label_file_name));
            
            %Create a label map for the new image based on the bilateral
            %filterred image
            idx = find(fim_1024bf~=-3);
            %Compute the confidence map
            newlblimbf = zeros(size(fim_1024bf));
            thresh =-0.025;
            candidateidx = find(lblim~=0);
            candidatelabel = (-1)*ones(size(candidateidx));
            %Smooth out fim to avoid sudden changes in the image
            %response
            glandidx = find(fim_1024bf(candidateidx)>thresh);
            candidatelabel(glandidx)=1; %Assign gland pixels
            newlblimbf(candidateidx)=candidatelabel;
            glandsize_file_name = strcat(label_name,'_morpbf.mat');
            lumenidx = find(newlblimbf==0); %Assign lumen pixels
            newlblimbf(lumenidx)=-1;
            %Clean the gland map a little bit
            newlblimbf = im2bw(newlblimbf,0);
            %Remove very small pores...
            newlbliminv = 1-newlblimbf;
            minporesize = 5;
            newlbliminv = bwareaopen(newlbliminv,minporesize,8); %Remove very small pores less...
            glandnoholes = 1-newlbliminv;
            %Slightly connect regions that are separated
            se = strel('disk',1);
            glandnoholes = imclose(glandnoholes,se);
            glandnoholes = bwareaopen(glandnoholes,5,8);
            newlblimbf = cast(glandnoholes,'single')*2.0-1.0;
            newlblimbf(lumenidx)=0;
            oldlblimbf = newlblimbf;
            [newlblimbf,gim,sim] = procseg(oldlblimbf); %Another processing step to remove all
            disp(['Calculating glandular dimensions of ' num2str(label_name)]);
                       
            
           if (~exist(glandsize_file_name,'file'))                
                %Caculate parameters for glands...
                glandmap = (newlblimbf==1);
                glandmap = imfill(glandmap,'holes');
                %Remove too small regions
                glandmap = bwareaopen(glandmap,minsizeinpixel);

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
                feat = [glandareamean,glandareamode,glandareamed,glandareastd,...
                glandpermean,glandpermode,glandpermed,glandperstd,...
                glandsoldmean,glandsoldmode,glandsoldmed,glandsoldstd,...
                glandstddistmean,glandstddistmode,glandstddistmed,glandstddiststd,...
                glandmaxdistmean,glandmaxdistmode,glandmaxdistmed,glandmaxdiststd,...
                glandcompactmean,glandcompactmode,glandcompactmed,glandcompactstd,...
                lumenareamean,lumenareamode,lumenareamed,lumenareastd,...
                lumenglandmean,lumenglandmode,lumenglandmed,lumenglandstd...
                ];
                figure(1);
                subplot(221);
                imagesc(newlblimbf);
                subplot(222);
                imagesc(glandmap);
                subplot(223);
                imagesc(lumenmap);
                subplot(224);
                plot(feat(end,:));
                drawnow;
                hold on;
               
                
               
                %Save the feature vector....
                save(strcat(morppath,glandsize_file_name),'feat','classidx');              
                nsamplesaved =  nsamplesaved+1;
                   
           else
                load(strcat(morppath,glandsize_file_name),'feat');
               
                    
           end
           dataclass(end+1,:)= classidx;
           totalfeat(end+1,:)=feat;
       end
    end
    close(h);
    save(strcat(morppath,'glandmorpbf.mat'),'totalfeat','dataclass');
end

function [newim,glandim,stromaimnew] = procseg(lblim)
%Post process the segmentation results
    %First, smooth out the gland image...
    lumenidx = find(lblim==0);
    lumenim = zeros(size(lblim));
    lumenim(lumenidx)=1;
    glandim = zeros(size(lblim));
    stromaim  = zeros(size(lblim));
    glandidx = find(lblim==1);
    stromaidx = find(lblim==-1);
    glandim(glandidx)=1;
    stromaim(stromaidx)=1;
    %glandim = bwareaopen(glandim, 16, 8); %Remove regions with too few pixels
    glandim = im2bw(glandim,0.5);     

    %Now, process the image of stroma
    stromaim = im2bw(stromaim,0.5);
    se = strel('disk',4); %Make sure the connection between two area is at least 8 pixels. We assume that all stroma is detected with high confidence
    stromaimnew = imclose(stromaim,se);%Further remove the holes and erode unecessary parts
    %Remove small object from the image...Edge information is not lost by
    %the closing operation....
    %through this operation
    minholesize = 50; %Remove isolated stroma area
    stromaimnew2 = bwareaopen(stromaimnew, minholesize, 4); %Remove regions with too few pixels
    %Remove gland regions that are two smalls
    stromaimnewinv = 1-stromaimnew2;
    minglandsize = 2000;
    stromaimnewinv2 = bwareaopen(stromaimnewinv, minglandsize, 4); %Remove regions with too few pixels
    stromaimnew3 = 1-stromaimnewinv2;
    
    %Remove too small lumen region
    lumenim1 = bwareaopen(lumenim,30);
    lumenidx = find(lumenim1==1);
    glandidx = setdiff(find(stromaimnew3==0),lumenidx);
    stromaidx = find(stromaimnew3==1);
    newim = zeros(size(lblim));
    newim(glandidx)=1;
    newim(stromaidx)=-1;   
end