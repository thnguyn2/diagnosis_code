function [featurename]=updateGlandSize(filenames,labelpath,diagpath,textonpath,mode,stromabandrad)
    %Last updated on April 24th
    %Added: g and histogram of texture for stroma description
    %This function computes the size of different glands in the biopsies
    %core...The glandsize as well as the distribution of glandsize will be
    %plotted...This code works on the 2048 x 2048 image. Otherwise, the
    %image has to be rescaled...
    %Compute the parameters of the biopsies using the groundtruth label...
    %To describe the stroma, we use the histogram of texton.
    %Inputs:
    %...
    %   stromabandrad: radius of the stroma band around the gland
    addpath(labelpath);
    addpath(diagpath);
    minsizeinpixel = 2000; %This is the minimum size for the gland area
    ntextons=50;
    h = waitbar(0,'Calculating size of the glands progress...');
    nfeatures = 25+2*ntextons;
    nsamplesaved = 0;
    feat = zeros(0,nfeatures);
    dataclass = zeros(0,1);
    ngrades = size(filenames,1);
    color_arr='rgbyck';
    retrain=0;
    featurename = cell(0,1);
   if ((~exist(strcat(diagpath,strcat('glandmorp_',mode,'_',num2str(stromabandrad),'.mat')),'file'))||(retrain==1))     
        
        for classidx=1:ngrades %Go through different classes
           nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
           for sampleIdx=1:nsamples
                waitbar(sampleIdx/nsamples,h,'Progress...')
                cur_file_name = filenames{classidx,1}{sampleIdx,1};
                %Check to see if the label file exist
                dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                slash_pos = strfind(cur_file_name,'\');
                label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                label_file_name = strcat(label_name,'_resized.mat');
                tx_file_name = strcat(label_name,'_texton.mat'); %Name for the texton map file

                disp(['Labeling: ' label_name ' ...']);
                smallim_name = strcat(cur_file_name(1:end-4),'_small.tif');
                phaseim = imread(smallim_name);
                phaseim = cast(phaseim,'single');
                load(strcat(textonpath,tx_file_name));
                if (strcmp(mode,'gt')==1)
                    load(strcat(labelpath,label_file_name));
                else
                    lblim = imread(strcat(cur_file_name(1:end-4),'_seg_multi_res_3072.tif'));                
                end
                [gf,th]= compute_gmod_features(lblim,phaseim,new_text_map,ntextons,stromabandrad); %Compute the map for gf
                [glandh] = compute_gland_features(lblim,phaseim,new_text_map,ntextons);

                     figure(3);
                     imagesc(lblim);
                     drawnow;
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


                    %Save the feature vector and the label...
                    feat(end+1,:) = [glandareamean,glandareamode,glandareamed,glandareastd,...
                    glandpermean,glandpermode,glandpermed,glandperstd,...
                    glandsoldmean,glandsoldmode,glandsoldmed,glandsoldstd,...
                    glandstddistmean,glandstddistmode,glandstddistmed,glandstddiststd,...
                    glandmaxdistmean,glandmaxdistmode,glandmaxdistmed,glandmaxdiststd,...
                    glandcompactmean,glandcompactmode,glandcompactmed,glandcompactstd,...
                    gf,th(:)',glandh(:)'];
                    dataclass(end+1,:)= classidx;
                    nsamplesaved =  nsamplesaved+1;
                    %figure(3);
                    %hold on;
                    clear areaarr;
                    clear boundaryarr;
                    clear centroidarr;
                    clear curbound2d;
                    clear glandcompactness;
                    clear glandh;
                    clear glandmap;
                    clear glandstddistancearr;
                    clear lblim;
                    clear new_text_map;
                    clear newim;
                    clear phaseim;
                    
                    
                    %plot(log10(feat(end,:)),color_arr(mod(classidx,6)+1));
                    %imagesc(feat);
            end
        end
        save(strcat(diagpath,strcat('glandmorp_',mode,'_',num2str(stromabandrad),'.mat')),'feat','dataclass');

   else
        load(strcat(diagpath,strcat('glandmorp_',mode,'_',num2str(stromabandrad),'.mat')),'feat','dataclass');
    end
   
    close(h);
 
end

function [gf,th] = compute_gmod_features(lblim,phaseim,textonim,ntextons,stromabandrad)
    %Return g_modified to the first argument and histogram of texton to the
    %2nd argument
    %This function compute the feature g_mod = 2*k0^2(1-g)=
    %<grad^2>/r/<var>_r
    se=strel('disk',stromabandrad);
    phaseim = phase_im_convert(phaseim/65536);
    %First, extract the stroma region
    stromamap = (lblim==2);
    stromamapinv = 1-stromamap;
    stromamapinv = bwareaopen(stromamapinv,1000);%Remove too small areas
    stromamap = 1-stromamapinv;    
    %To get the border, erode and then subtract
    stromamaperode=imerode(stromamap,se);
    stromamap = stromamap-stromamaperode;
  
    %Next, find the boundary of the stroma, get rid of the regions that is
    %the edge of the core
    newstromamap = (lblim==2);
    newstromamap = imfill(newstromamap,'holes');
    newmap = bwconvhull( newstromamap,'objects',4);
    enewstromamap = imerode(newmap,strel('disk',50));
    stromamap = stromamap.*enewstromamap;
    
    %Next, just use the regions around glands by certain radius
    glandmap = (lblim==1);
    glandmapdilate = imdilate(glandmap,strel('disk',stromabandrad));
    stromamap = stromamap.*glandmapdilate;
    figure(1);
    imagesc(stromamap);
    
    %Compute the variance of the phase
    phasevar = var(phaseim(stromamap==1));
    %Compute the gradient magnitude
    [gx,gy]=gradient(phaseim);
    gradmag = gx.^2+gy.^2;
    gradmag = mean(gradmag(stromamap==1));
    gf = gradmag/phasevar^2
    th=hist(textonim(stromamap==1),[1:ntextons])/sum(stromamap(:));
%     figure(2);
%     plot(th);
%     drawnow;
%     hold on;
    
end

function [th] = compute_gland_features(lblim,phaseim,textonim,ntextons)
    %th is the histogram of texton feature for the gland area
    glandmap = (lblim==1);
    th=hist(textonim(lblim==1),[1:ntextons])/sum(glandmap(:));
    %figure(2);
    %plot(th);
    %drawnow;
    %hold on;
    title('Histogram of texton for the gland regions');
    
end

function imrad = phase_im_convert(imdouble)
%This function perform a linear conversion where a double image (0 to 1 is converted into a phase image from -0.5*pi to 3.5*pi)
%Input: imdouble: an inpute image with grayscale values range from 0 to 1
%imrad: output image with phase values from -0.5 to 3.5
    min_val = -0.5;
    max_val = 3.5;
    imrad = (imdouble)*(max_val-min_val)+min_val;
end