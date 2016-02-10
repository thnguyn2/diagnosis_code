function diagnosis_feature_analysis
    clc;
    clear all;
    close all;
    addpath(strcat(cd(cd('..')),'\support'));
    datapath = 'G:\TMA_cores_and_diagnosis\';
    
    diagpath = 'G:\TMA_cores_and_diagnosis\diag\';
    featpath = 'G:\TMA_cores_and_diagnosis\feat\';%This is the folder for analysing diagnosis features
    disp('Finding the list of all labels');
    [filenames,glandnames,gradelist]=findFileNameFromROIs(datapath);
    save('filenameinfo.mat','filenames','glandnames','gradelist');
    mode = 'gt'; %Go with the segmentation based on the groundtruth
    %[textonfeat,classvect,samplename]=update_gland_texton_diagnosis_feature(filenames,strcat(datapath,'label\'),strcat(datapath,'texdir\'),mode,gradelist)
    %save('glandtextonfeat.mat','textonfeat','classvect','samplename');
    %[glandmorpfeat,classvect,samplename]=mean_size_analysis(filenames,strcat(datapath,'label\'),mode,gradelist);
    %save('glandmorpfeat.mat','glandmorpfeat','classvect','samplename');
    
    %[glandmorpfeat,classvect,samplename,indivglandfeat]=gland_morp_cancer_analysis(filenames,strcat(datapath,'label\'),mode,gradelist);%Compare the feature for the cancer group only
    %save('glandmorpfeat.mat','glandmorpfeat','classvect','samplename','indivglandfeat');
    load glandmorpfeat.mat;
    
     
%     %Show the histogram of features based on the Gleason grades     
    [mean2,mean3,mean4,mean5,std2,std3,std4,std5,cancergrade]=compute_feature_statistics_cancer_class(glandmorpfeat,classvect,true,samplename);%Draw the histogram of features for the cancer group
%     
%%     %Normalize w.r.t all classes so that the maximum value is 1
    mean_arr = [mean2;mean3;mean4;mean5];
    std_arr = [std2;std3;std4;std5];
    figure(5);
    barwitherr(std_arr', mean_arr');
    legend('G2','G3','G4','G5');
    set(gca,'XTickLabel',{'Gland arae', 'Peri ratio', 'GEcc','Non-Sold', 'Non-sold max', 'Stromainv','Mean welldef', 'Med welldef', 'Max welldef', 'h1','h2','h3','h4','h5'});
                      drawnow;

    %Do classification on 3 different class based on Gleason grades
    classify_cancerous_case(glandmorpfeat,cancergrade);
% %%==========================================================================
% 
% %     %Show the feature performance on Gleason's score
%     [meanlt7,meaneq7,meangt7,stdlt7,stdeq7,stdgt7,scorelist]=compute_feature_statistics_cancer_class_on_score(glandmorpfeat,classvect);%Draw the histogram of features for the cancer group
%     mean_arr = [meanlt7;meaneq7;meangt7];
%     std_arr = [stdlt7;stdeq7;stdgt7];
%     figure(4);
%     barwitherr(std_arr', mean_arr');
%     legend('LT7','EQ7','GT7');
%     set(gca,'XTickLabel',{'Area inhomo.', 'Non-solidity', 'Peri/cvx peri', 'Cribriform score'});
%     drawnow;
%     classify_cancerous_case_on_score(glandmorpfeat,scorelist);
% %=========================================================================
% 
% %     %Show the feature performance on for cancer/noncancer binary problem
%     [meanNM,meanHG,meanBP,meanCC,stdNM,stdHG,stdBP,stdCC,samplecond]=compute_feature_statistics_cancer_benign(glandmorpfeat,classvect);
%      mean_arr = [meanNM;meanHG;meanBP;meanCC];
%      std_arr = [stdNM;stdHG;stdBP;stdCC];
%      figure(4);
%      barwitherr(std_arr', mean_arr');
%      legend('NM','HG','BP','CC');
%      set(gca,'XTickLabel',{'Area inhomo.', 'Non-solidity', 'Peri/cvx peri', 'Cribriform score'});
% %=========================================================================
end

function [ls,g] = compute_ls_g(phasemap,wr)
    %Inputs:
    %   phasemap: a phase map for the values of the phase
    %   wr: radius of the filtering kernel that we use to compute the ls
    %   and sg 
    %Outputs:
    %   ls, g; scattering mean freepath and the anisotropy factor
    %Formulae:
    %   ls = 1/<phase_variance over the windows>
    %   sg = 1 - <gradient magnitute>/<phase variance>
    %This function compute the ls and g for the phase image
    %Compute the ls map: scattering mean freepath
    phasemap = cast(phasemap,'single');
    maxval = 65535;
    minval = 0;
    h = fspecial('disk',wr);
    %First, compute the mean over the wirndows
    meanphase = imfilter(phasemap,h,'same');
    errorphase = phasemap - meanphase;%Compute the difference betweeen each pixel and the surrounding mean
    errorphase = errorphase.^2;%Square the error
    varphase = imfilter(errorphase,h,'same');%Compute the final average and invert it
    ls = 1./varphase;
    
    %Compute the gradient map
    [gx,gy]=gradient(phasemap);
    gmag2 = gx.^2+gy.^2; %Gradient magnitute squared
    %Now, filter the gradient magnitude
    gmag2avg = imfilter(gmag2,h,'same');
    g = 1-gmag2avg./varphase; %This is a measure of anisotropy
  
end
function [feat,dataclass,samplename,indivglandfeat]=gland_morp_cancer_analysis(filenames,labelpath,mode,gradelist)
%Analize the morphological feature of the glands for all the cancer cases
%Outputs:
%   indivglandfeat: a vector of cells where each cell is a 2D array of Ni x
%   P where Ni is the number glands and P is the number of features
      dispen = 0;
      P = 17; %Number of feature for diagnosis on the whole core....
      feat = zeros(0,P);
      indivglandfeat = cell(0,1);
      dataclass = cell(0,1); %This is the class of the sample.
      ngrades = size(filenames,1);
      samplename = cell(0,1);
      minglandsizeinpixel = 10000; %Define the minimum size of the gland in pixels. This is to avoid concerning about to small glands
      maxlumenbandwidthsurroundinggland=15; %This is the maximum width of the lumen band surrounding each gland. If there is lumen, it will reduce the stroma probability, hence increases the well definess of gland.
      boundarywidth=20;             %This is the width of the boundary for handling gland's boundary
      glandboundarywidth = 20;
      glandclearthreshold = 1; %If a segment that needs to be added to a gland is longer than this ratio times the perimeter of the gland's convexhull, the gland will not be used in the calculation. Set this factor to 1 to include all glands
      addedlumeneccenthresh = 0.9; %To-be-added lumen areas over this value will not be added
      minaddedlumenarea = 10000; %Area over this area will be added irrespective of there eccentricity
      minlumensizeinpixel = 5000;
      minnonlumensize = 100000;   %A threshold to get rid of all the details (debris)
      minglandboundaryifcut = 20; %Thin width of the glands....
      mindistanceforglandjoint = 40;
      minstromawidth = 10;
      lsgwr = 20;
      for classidx=1:ngrades %Go through different classes
           if (~ismember(gradelist(classidx),{'NC','SK'})) %Only check the cancer case
               disp(['Working on:' gradelist(classidx)]);
               nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
               for sampleIdx=1:min(nsamples,nsamples)
               %for sampleIdx=4:4
                 
                    
                    cur_file_name = filenames{classidx,1}{sampleIdx,1};
                    %Check to see if the label file exist
                    dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                    slash_pos = strfind(cur_file_name,'\');
                    label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                    label_file_name = strcat(label_name,'_resized.mat');
                    disp(['Working on : ' label_name ' ...']);
                    phaseim = imread(strcat(cur_file_name(1:end-4),'_small.tif'));%This is the phase map to compute ls and g
                    %Read the probability map of stroma
                    pmapim = imread(strcat(cur_file_name(1:end-4),'_seg_ws90_fuzzy.tif'));
                    if (strcmp(mode,'gt')==1)
                        load(strcat(labelpath,label_file_name));
                       
                    else
                        lblim = imread(strcat(cur_file_name(1:end-4),'_seg_multi_res_3072.tif'));                
                    end
                    %Prepare the ls and g map for the features
                    [ls,g] = compute_ls_g(phaseim,lsgwr);
                    
                    %Generate a good map of all the label. Handle the
                    %following cases:
                    % +Debris.
                    % +Very small glands => Eliminate and set it to stroma.
                    % ***The small glands must be well-formed***Otherwise,
                    % just eliminate them.
                    % +Boundary glands that may be over estimated for the features...
                    % +Very small regions inside the gland
                    % +A small sheet of epithelial cell broken that leads
                    % to wrong estimation of gland's shape
                    % + Lumen identification: A gland very thin and lumen invade it.
                    % +Gland's boundary specification
                    %The first case can be done by finding all the glands. Circumscribe them with
                    %a convex hull. Find the interesection for the convex
                    %hull of each gland with the non-lumen map. If there
                    %is a section of the gland that exist in the non-lumen
                    %map, assign that area to lumen and circumscribe it
                    %with a gland boundary. 
                    
                  
                    
                    %Zero processing- Debris rejection
                    nonlumenmap = lblim~=0;
                    nonlumenmap = imopen(nonlumenmap,strel('disk',5)); %Remove very small connection of debris
                    areaarr = regionprops(nonlumenmap,'Area');
                    nregions = size(areaarr,1);
                    nonlumenidxlist = regionprops(nonlumenmap,'PixelIdxList');
                    for regionidx=1:nregions
                        if (areaarr(regionidx).Area<minnonlumensize)
                            nonlumenmap(nonlumenidxlist(regionidx).PixelIdxList)=0;%Clear very small area
                            lblim(nonlumenidxlist(regionidx).PixelIdxList)=0;   %Update the label map
                        end
                    end
                    
                    
                    %First processing, get rid of very small glands
                    %Expand the stroma to cut out small connection between
                    %glands
                    stromamap = lblim==2;
                    stromamap = imclose(stromamap,strel('disk',minstromawidth));
                    stromapixidx=find(stromamap==1);
                    lblim(stromapixidx)=2;
                    glandmap = lblim==1;
                    areaarr = regionprops(glandmap,'Area');
                    nregions = size(areaarr,1);
                    glandpixellist = regionprops(glandmap,'PixelIdxList');
                    for regionidx=1:nregions
                        if (areaarr(regionidx).Area<minglandsizeinpixel)
                            glandmap(glandpixellist(regionidx).PixelIdxList)=0;%Clear very small area
                            lblim(glandpixellist(regionidx).PixelIdxList)=2;   %Update the label map
                        end
                    end
                    
                    %Second processing - Boundary glands   
                    nonlumenmap = imfill(nonlumenmap,'holes');
                    nonlumenidx = find(nonlumenmap==1);
                    
                    %Find the core's boundary                    
                    nonlumenmap2 = imerode(nonlumenmap,strel('disk',boundarywidth));
                    boundarymap = nonlumenmap-nonlumenmap2;
                    boundarypixidx = find(boundarymap==1);
                    conveximages = regionprops(glandmap,'ConvexImage'); %3
                    boundingboxes  = regionprops(glandmap,'BoundingBox');
                    glandpixellist = regionprops(glandmap,'PixelIdxList');   
                    glandmap = cast(glandmap,'uint8');
                    nregions = size(conveximages,1);
                    nrows=size(lblim,1);
                    addedlumenmap = zeros(size(lblim));
                    
                    for regionidx=1:nregions
                        if (~isempty(intersect(boundarypixidx,glandpixellist(regionidx).PixelIdxList)))
                                r1 = round(boundingboxes(regionidx,1).BoundingBox(2));%Get the coordinates of the current glands
                                c1 = round(boundingboxes(regionidx,1).BoundingBox(1));
                                [glandrowidx,glandcolidx]=find(conveximages(regionidx).ConvexImage); %Find the index of the pixels in the circumscribe
                                cvxglandrowidx = glandrowidx+r1-1;
                                cvxglandcolidx = glandcolidx+c1-1;
                                cvxglandidx = cvxglandrowidx+(cvxglandcolidx-1)*nrows;
                                outsidelumenidx = setdiff(cvxglandidx,nonlumenidx);%Set all the convex area outside the core to
                                addedlumenmap(outsidelumenidx)=1;
                                %Another steps, find the added length of a
                                %segement to close the gland
                                tempcvxim = conveximages(regionidx).ConvexImage;
                                [glandboundaryrowidx,glandboundarycolidx] = find((tempcvxim-imerode(tempcvxim,strel('disk',glandboundarywidth)))==1);%Find the coordinates of the boundary for the convex area surrounding the gland
                                glandboundaryrowidx = glandboundaryrowidx + r1-1;
                                glandboundarycolidx = glandboundarycolidx + c1-1;
                                glandboundaryidx = glandboundaryrowidx+( glandboundarycolidx-1)*nrows;
                                addedlumensegidx = setdiff(glandboundaryidx,nonlumenidx);
                                curaddedsegratio = length(addedlumensegidx)/length(glandboundaryidx);
                                if (curaddedsegratio>glandclearthreshold)
                                    addedlumenmap(outsidelumenidx)=0; %Do not added the lumen map
                                    glandmap(glandpixellist(regionidx).PixelIdxList)=0; %Clear the gland map
                                    lblim(glandpixellist(regionidx).PixelIdxList)=2; %Clear the label map
                                end
                               
                        end
                    end
           
                        addedlumenmap = bwareaopen(addedlumenmap, minlumensizeinpixel);%Avoid adding lumen area that area too small
                        %Avoiding adding lumen areas that is too thin...
                        addedlumeneccent = regionprops(addedlumenmap,'Eccentricity');
                        addedpixelidx = regionprops(addedlumenmap,'PixelIdxList');  
                        addedlumenarea = regionprops(addedlumenmap,'Area');
                        naddedregions = length(addedlumeneccent);
                        for regionidx=1:naddedregions
                            if ((addedlumeneccent(regionidx).Eccentricity>addedlumeneccenthresh)&(addedlumenarea(regionidx).Area<minaddedlumenarea))
                                addedlumenmap(addedpixelidx(regionidx).PixelIdxList)=0; %Clear the added region
                            end
                        end
                        %...............................................
                        addedglandboundary = addedlumenmap-imerode(addedlumenmap,strel('disk',5));
                        %Update the glandmap and the label map
                        addglandboundaryidx = find(addedglandboundary==1);
                        glandmap(addglandboundaryidx)=1;
                        lblim(addglandboundaryidx)=1;
                    
                    
                    %Thrid processing - removing small luminal area inside
                    %the gland
                    glandmapfilled = imfill(glandmap,'holes');
                    invglandmap = cast((1-glandmap),'uint8').*cast((glandmapfilled==1),'uint8');
                    invglandmap = bwareaopen(invglandmap,minlumensizeinpixel);
                    glandmap = cast((1-invglandmap),'uint8').*cast((glandmapfilled==1),'uint8');
                    glandpixidx = find(glandmap==1);
                    lblim(glandpixidx)=1;
                    
                    %Forth processing -  Handling very thing glands
                    %boundary. This is for specifying the lumen in fifth
                    %step. Once the lumen is specified, we repeat step 4 to
                    %specify gland's boundary
                    
                    glandmappreserved = glandmap;%Preserve the gland map and the label map
                    lblimpreserved = lblim;
                    
                    glandmap = glandmap>=1;%Convert back to the binary case
                    glandidxlist = regionprops(glandmap,'PixelIdxList');
                    nregions = size(glandidxlist,1);
                    for regionidx=1:nregions
                        curglandidx1 = glandidxlist(regionidx).PixelIdxList;              %Find all 2D coordinates of the current gland
                        [curglandrowidx,curglandcolidx]=ind2sub(size(lblim),curglandidx1);%of the current gland
                        r1 = min(curglandrowidx);
                        c1 = min(curglandcolidx);
                        r2 = max(curglandrowidx);
                        c2 = max(curglandcolidx);
                        bbnrows = r2-r1+1;
                        bbncols = c2-c1+1;
                        tempglandim = zeros(bbnrows,bbncols);
                        curglandrowidx = curglandrowidx+1-r1;
                        curglandcolidx = curglandcolidx+1-c1;
                        curglandidx = curglandrowidx + (curglandcolidx-1)*bbnrows;
                        tempglandim(curglandidx)=1;
                        tempglandim = imclose(tempglandim,strel('disk', minglandboundaryifcut));%Close very small boundary, act on a single gland and leave other glands alone
                        tempglandim = imfill(tempglandim,'holes'); %We don't want to close anything inside the gland
                        tempglandnbound = tempglandim - imerode(tempglandim,strel('disk',3));%Just care about the boundary since we don't want
                        %to squeeze the lumen area inside
                        [curglandrowboundidx,curglandcolboundidx]=find(tempglandnbound==1);
                        curglandrowboundidx = curglandrowboundidx+r1-1;%Convert to global coordinates
                        curglandcolboundidx = curglandcolboundidx+c1-1;%Convert to global coordinates
                        curglandidx2 = curglandrowboundidx+(curglandcolboundidx-1)*nrows; %Coordinates of the new gland
                        addedglandidx=setdiff(curglandidx2,curglandidx1);
                        glandmap(addedglandidx)=1;
                        lblim(addedglandidx)=1;
                    end
                    %Fifth processing step - Lumen identification
                    glandidxlist = regionprops(glandmap,'PixelIdxList');
                    nregions = size(glandidxlist,1);
                    for regionidx = 1:nregions
                         curglandidx1 = glandidxlist(regionidx).PixelIdxList;              %Find all 2D coordinates of the current gland
                        [curglandrowidx,curglandcolidx]=ind2sub(size(lblim),curglandidx1);%of the current gland
                        r1 = min(curglandrowidx);
                        c1 = min(curglandcolidx);
                        r2 = max(curglandrowidx);
                        c2 = max(curglandcolidx);
                        bbnrows = r2-r1+1;
                        bbncols = c2-c1+1;
                        tempglandim = zeros(bbnrows,bbncols);
                        curglandrowidx = curglandrowidx+1-r1;
                        curglandcolidx = curglandcolidx+1-c1;
                        curglandidx = curglandrowidx + (curglandcolidx-1)*bbnrows;
                        tempglandim(curglandidx)=1;
                        tempglandimfill = imfill(tempglandim,'holes');
                        roilblim =  lblim(:,c1:c2);
                        roilblim=roilblim(r1:r2,:);
                        %Specify the lumen area
                        lumenmap = roilblim==0;
                        lumenmap = lumenmap.*tempglandimfill;
                        %Eliminate very small luminal regions
                        lumenmap = bwareaopen(lumenmap,minlumensizeinpixel);
                        lumenmap = imfill(lumenmap,'holes');%Make sure there is no debris inside the lumen
                        [lumenrowrelidx,lumencolrelidx]=find(lumenmap==1);
                        lumenpixidx = (lumencolrelidx+c1-2)*nrows+lumenrowrelidx+r1-1;
                        lblim(lumenpixidx)=-1;
                    end
                    
                    %Now, update found lumen coordinates to the gland map
                    %Sixth step: specifying the gland's boundary
                    lumenidx = find(lblim==-1);
                    lblimpreservedglandfilled=lblimpreserved;
                    lblimpreservedglandfilled(lumenidx)=1;%This is the final map of the glands.
                    glandmap = lblimpreservedglandfilled==1;
                    glandmapboundary = glandmap-imerode(glandmap,strel('disk',5));%Final map of the gland's boundary
                    %Return the correct lumen ara
                    lblimpreserved(lumenidx)=-1;
                    glandboundarypixidx = find(glandmapboundary==1);
                    lblimpreserved(glandboundarypixidx)=1;
                    lblim = lblimpreserved;
                    glandmap = lblim==1;
                    noholeglandmap = imfill(glandmap,'holes');
                    
                    if (dispen)

                        figure(1);
                        imagesc(glandmap);drawnow;title('Gland map')
                        figure(2);
                        imagesc(lblim);drawnow;title('label map'); 

                        %Generate a 6 shade color map


                        figure(3);
                        subplot(121);
                        imagesc(pmapim);drawnow;title('Probability map'); colormap bone;colorbar
                        subplot(122);
                        imagesc(histeq(phaseim));
                        colorbar;colormap gray;
                        title(label_name);
                        drawnow;

                    end
                    %Step 2: compute different features
                    %   +Specify the centroid of all the cores (for
                    %   display)
                    %   +Use the gland size features
                    %   +Gland's perimeters/perimeter of the fitting ellipse
                    %   +Gland's eccentricity
                    %   +Gland's solidity
                    %   +Invasion score
                    %   +Well-definess score.
                    %   +Lumen's area/gland's area
                    %   +Cribriform score
                    %   +Anisotropy map...
                    %   +Smoothness of gland's boundary
                    
                    %Care about the ls and g map for the gland only
                    ls = ls.*glandmap;
                    g = g.*glandmap;
%                     figure(6);
%                     subplot(121); imagesc(ls);title('LS');colorbar
%                     subplot(122); imagesc(g);title('SG');colorbar
%                     
    
    
                    disp('Computing parameters of the cores...');
                    %0th step: get the coordinates of all the glands for
                    %displaying the features' value and gland size
                    centroidarr = regionprops(glandmap,'Centroid');
                    nregions=size(centroidarr,1);
                    xcoord = zeros(nregions,1);
                    ycoord = zeros(nregions,1);
                    for regionidx=1:nregions
                        cur_cent = (centroidarr(regionidx,1).Centroid);    %First column is X, second column is Y 
                        xcoord(regionidx,1)=cur_cent(1);
                        ycoord(regionidx,1)=cur_cent(2);                        
                    end
                    
                    
                    curF=zeros(nregions,6);%Add the new feature to the feature map
                   
                    %1st feature: glandsize features
                    glandsizearr = regionprops(noholeglandmap,'Area');
                    gsizefeat = zeros(nregions,1);
                    for regionidx=1:nregions
                        gsizefeat(regionidx,1)=glandsizearr(regionidx,1).Area;
                    end
                    gsizemean = 1-mean(gsizefeat)/max(gsizefeat);
                    curF(:,1)=gsizefeat/max(gsizefeat);
                                       
                    %2nd feature: gland's perimeter/perimetter of the
                    %fitting ellipse
                    periarr = regionprops(glandmap,'Perimeter'); %3
                    majaxislength = regionprops(glandmap,'MajorAxisLength');
                    minaxislength = regionprops(glandmap,'MinorAxisLength');
                    nregions = size(periarr,1);%Get the number of regions
                    perifeat= zeros(nregions,1);
                    for regionidx=1:nregions
                        perifeat(regionidx,1)=periarr(regionidx).Perimeter/(pi*sqrt(majaxislength(regionidx).MajorAxisLength^2+...
                            minaxislength(regionidx).MinorAxisLength^2));
                    end
                    perimean = mean(perifeat.*gsizefeat)/mean(gsizefeat);
                    curF(:,2)=perifeat;
                    
                    %3rd feature: gland's eccentricity
                    geccarr = regionprops(noholeglandmap,'Eccentricity'); %Get the eccentricity of the lumen area
                    geccfeat = zeros(nregions,1);
                    for regionidx=1:nregions
                        geccfeat(regionidx)=geccarr(regionidx).Eccentricity;
                    end
                    %geccmean = mean(geccfeat.*gsizefeat)/mean(gsizefeat);
                    geccmean = mean(geccfeat);
                    
                    curF(:,3)=geccfeat;
                    
                    %4th feature: gland's solidity
                    soldarr = regionprops(noholeglandmap,'Solidity');
                    soldfeat = zeros(nregions,1);
                    for regionidx=1:nregions
                        soldfeat(regionidx)=soldarr(regionidx).Solidity;
                    end
                    soldfeat = 1-soldfeat; %We want to look at the non-solidity of the gland
                    %soldmean = mean(soldfeat.*gsizefeat)/mean(gsizefeat);
                    soldmean = mean(soldfeat);
                    soldmax = max(soldfeat);
                    curF(:,4)=soldfeat;
                    
                    %5th feature: stromainvasion score
                    sivsfeat = zeros(nregions,1);
                    conveximages = regionprops(glandmap,'ConvexImage'); %3
                    boundingboxes  = regionprops(glandmap,'BoundingBox');
                    warning('off');
                    for regionidx=1:nregions
                               currentconvexim = conveximages(regionidx).ConvexImage;%Get a faction of the label map inside the convex hull
                               currentbox = boundingboxes(regionidx).BoundingBox;
                               r1 = currentbox(2);
                               c1 = currentbox(1);
                               nr = currentbox(4);
                               nc = currentbox(3);
                               r2 = r1+nr-1;
                               c2 = c1+nc-1;
                               smallconvexgland = currentconvexim.*lblim(r1:r2,c1:c2);
                               row_centr = ycoord(regionidx)-r1;%Relative coordinate of the centroid
                               col_centr = xcoord(regionidx)-c1;
                               radius =   minaxislength(regionidx).MinorAxisLength/2;%This is the radius for computing the stromal invasion
                               sivsfeat(regionidx)= compute_stroma_invasion(row_centr,col_centr,smallconvexgland,radius,sum(sum(currentconvexim)));%Compute the stromal invasion feature for each gland   
                    end
                    sivsmean = mean(sivsfeat.*gsizefeat)/mean(gsizefeat);
                    curF(:,5)=sivsfeat;
                    
                    
                    %6th feature: the probability on well-definess of the
                    %gland. it ts the average probability of the core.
                  %  glandmapwidthboundary = imclose(glandmap,strel('disk',maxlumenbandwidthsurroundinggland));
                    welldef = zeros(nregions,1);
                    glandidxlist = regionprops(glandmap,'PixelIdxList');
                    for regionidx=1:nregions
                        curglandidx1 = glandidxlist(regionidx).PixelIdxList;              %Find all 2D coordinates of the current gland
                        [curglandrowidx,curglandcolidx]=ind2sub(size(lblim),curglandidx1);%of the current gland
                        r1 = min(curglandrowidx);
                        c1 = min(curglandcolidx);
                        r2 = max(curglandrowidx);
                        c2 = max(curglandcolidx);
                        bbnrows = r2-r1+1;
                        bbncols = c2-c1+1;
                        tempglandim = zeros(bbnrows,bbncols);
                        curglandrowidx = curglandrowidx+1-r1;
                        curglandcolidx = curglandcolidx+1-c1;
                        curglandidx = curglandrowidx + (curglandcolidx-1)*bbnrows;
                        tempglandim(curglandidx)=1;
                        tempglandimwithlumenbound = imclose(tempglandim,strel('disk',maxlumenbandwidthsurroundinggland));
                        [glandrowidx,glandcolidx]=find(tempglandimwithlumenbound==1);
                        glandrowidx = glandrowidx + r1-1;
                        glandcolidx = glandcolidx + c1-1;
                        curglandidx2 = glandrowidx+(glandcolidx-1)*nrows; %Coordinates of the new gland
                                           
                        welldef(regionidx) = mean(pmapim(curglandidx2));
                    end
                    maxwelldef = max(welldef);
                    meanwelldefscore = mean(welldef);
                    medwelldef = median(welldef);
                    histbinval = linspace(0,0.5,8);
                    curhist=hist(welldef,histbinval)/length(welldef);
                    
                    
                    curF(:,6) =welldef;
                    
                    indivglandfeat{end+1,1}=curF;
                    dataclass(end+1,:)= gradelist(classidx);
                    samplename{end+1}=label_name;
                    feat(end+1,:)=[gsizemean perimean geccmean soldmean soldmax sivsmean meanwelldefscore medwelldef maxwelldef curhist(:)'];
                    warning('on');
                       
                    if (dispen)
                        %Display the features on the top of the image
                        for regionidx=1:nregions
                             figure(2);
                             text(xcoord(regionidx), ycoord(regionidx), sprintf('%1.3f',welldef(regionidx)), 'Color', 'r', ...
                             'FontWeight', 'bold');                         

                        end
                        disp(['Current feature: ' num2str(curhist(4))])
                        [mean2,mean3,mean4,mean5,std2,std3,std4,std5,cg,c2f,c3f,c4f,c5f]=compute_feature_statistics_cancer_class(feat,dataclass,false,samplename);

                        %Normalize w.r.t all classes so that the maximum value is 1
                         mean_arr = [mean2;mean3;mean4;mean5];
                         std_arr = [std2;std3;std4;std5];

                         figure(7);
                         barwitherr(std_arr', mean_arr');
                         legend('G2','G3','G4','G5');                    
                         drawnow;
                         hold on;
                         plot(feat(end,:),'-xr');
                         hold off;
                          set(gca,'XTickLabel',{'Gland arae', 'Peri ratio', 'GEcc','Non-Sold', 'Non-sold max', 'Stromainv','Mean welldef', 'Med welldef', 'Max welldef', 'h1','h2','h3','h4','h5'});
                          drawnow;

                          %Convert the feature to the histogram
                          featidx=13;
                          curfeat_arr = cell(4,1);
                          curfeat_arr{1} = c2f(:,featidx); 
                          curfeat_arr{2} = c3f(:,featidx); 
                          curfeat_arr{3} = c4f(:,featidx); 
                          curfeat_arr{4} = c5f(:,featidx); 

                          %Draw the scatter plot of features
                          n_arr = zeros(4,1);
                          for gradeidx=1:4
                               n_arr(gradeidx) = length(curfeat_arr{gradeidx});
                          end

                          figure(4);
                          colorarr = 'bgrm';
                          for gradeidx=1:4
                              if (n_arr(gradeidx)>0)
                                for sampleidx=1:n_arr(gradeidx)
                                    %plot(gradeidx,curfeat_arr{gradeidx}(sampleidx),strcat('o',colorarr(gradeidx)),'linewidth',2);
                                    figure(4);
                                    plot(curhist,strcat('-o',colorarr(gradeidx)),'linewidth',2);
                                    xlabel('Bin idx');
                                    grid on;
                                    hold on;
                                    drawnow;
                                    figure(5);
                                    [gsizefeatsorted,sortidx] = sort(gsizefeat);
                                    welldefsorted = welldef(sortidx);
                                    plot(gsizefeatsorted,welldefsorted,strcat('-o',colorarr(gradeidx)),'linewidth',2);
                                    xlabel('Gland size');
                                    ylabel('Welldef score');
                                    grid on;
                                    hold on;
                                    drawnow;
                                end
                              end
                          end
                     end
                      disp('Done');
    
               end
            end
   
      end
end
function [si] = compute_stroma_invasion(row_centr,col_centr,lblim,radius,convexglandarea)
%Compute the stromal invasion features
%row_centr,col_centr: row and column coordinates of the centroid
%lblim: a small label image for each gland
%convexglandarea: the area of the convexhull surrounding each gland
%si: stroma invasion index
  
    weightmap = zeros(size(lblim));
    [nrows,ncols]=size(lblim);
    %Compute the distance between each stromal pixels and the centroid
    [rowidx,colidx]=find(lblim==2);
    pixidx = (colidx-1)*nrows+rowidx;
    dist_sqr = (rowidx-row_centr).^2+(colidx-col_centr).^2;%Distance between every pixels to the centroid
    weight_vect =exp(-dist_sqr/2/radius.^2);
    weightmap(pixidx)=weight_vect;
    si = sum(weightmap(:))/convexglandarea;  
end
function classify_cancerous_case(glandmorpfeat,classvect)
    %This function perform cancer classification by forming several binary
    %classification and 1 multi-class classification
    %Classify between cancer grade 2 and 3
    ctype = 'ensemble';
    nround = 3;
    %glandmorpfeat=glandmorpfeat(:,[1:3 5:7 9:12]);
%  
     kval = 5;%-1 mean leave 1 out cross-validation
%     disp('Grade 2 vs grade 3 - confusion matrix:')
%     idx = union(find(classvect==2),find(classvect==3));
%     data = glandmorpfeat(idx,:);species = classvect(idx,:);
%     cf_matrix=classifywithcv(data,species,nround,kval,ctype);
%     disp(num2str(cf_matrix));
%     
%      
%     disp('Grade 2 vs grade 4 - confusion matrix:')
%     idx = union(find(classvect==2),find(classvect==4));
%     data = glandmorpfeat(idx,:);species = classvect(idx,:);
%     cf_matrix=classifywithcv(data,species,nround,kval,ctype);
%     disp(num2str(cf_matrix));
%     
%     disp('Grade 2 vs grade 5 - confusion matrix:')
%     idx = union(find(classvect==2),find(classvect==5));
%     data = glandmorpfeat(idx,:);species = classvect(idx,:);
%     cf_matrix=classifywithcv(data,species,nround,kval,ctype);
%     disp(num2str(cf_matrix));
    
    disp('Grade 3 vs grade 4 - confusion matrix:')
    idx = union(find(classvect==3),find(classvect==4));
    data = glandmorpfeat(idx,:);species = classvect(idx,:);
    cf_matrix=classifywithcv(data,species,nround,kval,ctype);
    disp(num2str(cf_matrix));
    
    disp('Grade 3 vs grade 5 - confusion matrix:')
    idx = union(find(classvect==3),find(classvect==5));
    data = glandmorpfeat(idx,:);species = classvect(idx,:);
    cf_matrix=classifywithcv(data,species,nround,kval,ctype);
    disp(num2str(cf_matrix));
    
    
%     disp('Grade 4 vs grade 5 - confusion matrix:')
%     idx = union(find(classvect==4),find(classvect==5));
%     data = glandmorpfeat(idx,:);species = classvect(idx,:);
%     cf_matrix=classifywithcv(data,species,nround,kval,ctype);
%     disp(num2str(cf_matrix));
   
 
    disp('Grade 2,3,4 - confusion matrix:')
    idx = [find(classvect==2);find(classvect==3);find(classvect==4)];
    data = glandmorpfeat(idx,:);species = classvect(idx,:);
    cf_matrix=classifywithcv(data,species,nround,kval,ctype);
    disp(num2str(cf_matrix));
    
    disp('Grade 2,3,4,,5 - confusion matrix:')
    idx = [find(classvect==2);find(classvect==3);find(classvect==4);find(classvect==5)];
    data = glandmorpfeat(idx,:);species = classvect(idx,:);
    cf_matrix=classifywithcv(data,species,nround,kval,ctype);
    disp(num2str(cf_matrix));
end%Classify cancer cases only
function classify_cancerous_case_on_score(glandmorpfeat,scorelist)
    %This function perform cancer classification by forming several binary
    %classification and 1 multi-class classification
    %Classify between cancer grade 2 and 3
    ctype = 'ensemble';
    nround = 10;
    glandmorpfeat=glandmorpfeat(:,1:3);
    kval = 5;%-1 mean leave 1 out cross-validation
    disp('Score lt7 vs score 7 - confusion matrix:')
    idx = union(find(scorelist==6),find(scorelist==7));
    data = glandmorpfeat(idx,:);species = scorelist(idx,:);
    cf_matrix=classifywithcv(data,species,nround,kval,ctype);
    disp(num2str(cf_matrix));
    
     
    disp('Score 7 vs score gt 7 - confusion matrix:')
    idx = union(find(scorelist==7),find(scorelist==8));
    data = glandmorpfeat(idx,:);species = scorelist(idx,:);
    cf_matrix=classifywithcv(data,species,nround,kval,ctype);
    disp(num2str(cf_matrix));
    
    disp('Score lt7 vs score gt7 - confusion matrix:')
    idx = union(find(scorelist==6),find(scorelist==8));
    data = glandmorpfeat(idx,:);species = scorelist(idx,:);
    cf_matrix=classifywithcv(data,species,nround,kval,ctype);
    disp(num2str(cf_matrix));

    disp('Score lt7, 7, gt7 - confusion matrix:')
    idx = [find(scorelist==6);find(scorelist==7);find(scorelist==8)];
    data = glandmorpfeat(idx,:);species = scorelist(idx,:);
    cf_matrix=classifywithcv(data,species,nround,kval,ctype);
    disp(num2str(cf_matrix));
end%Classify based on Gleason's score

function [meanval,stdval,skewval,kurtval,maxval]=get_pdf_moment(X,donorm,glandsizearr)
    %Compute the moments of the distribution for the variable X
    %donorm: if true: do normalization w.r.t the sizes of glands
    %glandsizearr: an array of gland size
    if (nargin==1)
        donorm=false;    
    end
    meanval = mean(X.*glandsizearr)/mean(glandsizearr);
    kurtval = kurtosis(X(:));
    stdval = std(X(:));
    skewval = skewness(X(:));
    maxval = max(X(:));
  
    if (donorm==false)
        meanval = mean(X);
        %Use the empirical distribution to weight the sample
    end
    
end
function [cftotalnew]=classifywithcv(data,species,nround,kval,ctype)
    %Perform nround of classifycation with k-fold cross-validataion with
    %classifier of type ctype. The number k is specified by kval.
    %ctype belongs to 'bayes'/'svm'/'rf'/'disc'/'nn'
    %data: feature matrix where each row is an observation, each column is
    %a feature
    nclass = length(unique(species));
    cftotal = zeros(nclass,nclass);
    nfeat=size(data,2);
    for roundidx=1:nround %Each round start with a different random seeds
         rng(roundidx*20); %Random seed generator
         if (kval~=-1)
             CVO = cvpartition(species,'kfold',kval); %Go with 4 since we have at least 1 sample for each class
             ntestset = CVO.NumTestSets;
         else
             ntestset = size(data,1);
         end
         err = zeros(ntestset,1);
         cp_arr=zeros(ntestset,1);
         cfmat = cell(ntestset,1);
        
         for i = 1: ntestset;
              if (kval~=-1)
                  trIdx = CVO.training(i);
                  teIdx = CVO.test(i);
              else
                  teIdx=zeros(size(species));
                  trIdx=ones(size(species));
                  teIdx(i)=1;
                  trIdx(i)=0;
              end
              switch ctype
                  case 'bayes'
                      %Bayesian classifier
                      c = NaiveBayes.fit(data(trIdx,:),species(trIdx,:),'Prior','uniform');%Use uniform prior for each class 
                      yout = c.predict(data(teIdx,:));
                  case 'svm'
                      %Svm classifier - works for data after PCA reduction
                      %only since the number of observation is smaller than
                      %the number of features in our problem
                      c=0.05; %C needs to be small to be general enough. Unbalance data should be handled automatically
                      rbf_sigma=0.8; %This constant needs to be large enough. If too small=> Nearly linear decision boundary. Too large=> overfitting
                      if (nfeat==2)
                          svmstruct = svmtrain(data(trIdx,:),species(trIdx,:),'kernel_function','rbf','rbf_sigma',rbf_sigma...
                              ,'boxconstraint',c,'showplot','true');
                      else
                          svmstruct = svmtrain(data(trIdx,:),species(trIdx,:),'kernel_function','rbf','rbf_sigma',rbf_sigma...
                              ,'boxconstraint',c);

                      end
                      yout = svmclassify(svmstruct,data(teIdx,:));
                      
                  case 'rf'
                      ntrees = 100;
                      rfstruct = RTrees;
                      rfstruct.train(data(trIdx,:),species(trIdx,:),'MaxNumOfTreesInTheForest',ntrees,...
                        'ForestAccuracy',0.05,'MaxDepth',8,'Prior',[1 1],'VarType','Categorical','CalcVarImportance',1); %Putting more weights on stroma.
                      yout = rfstruct.predict(data(teIdx,:));
                      evalfeatimp=0;
                      if (evalfeatimp)
                          %Evaluate the importance of features
                          varimp = rfstruct.getVarImportance();
                          accum_varimp = (accum_varimp*(i-1)+varimp)/i; %Accumulated feature importance
                          [U,V]=sort(accum_varimp,'descend');
                          feat1idx=V(1);feat2idx=V(2);
                          gscatter(data(:,feat1idx),data(:,feat2idx),cannormspecies,'rgbcy');
                          xlabel(num2str(feat1idx));
                          ylabel(num2str(feat2idx));
                          drawnow
                      end
                  case 'ensemble' %Go with the tree bagging classifier
                      t = templateTree('minleaf',2);%Create a deep tree.
                      B = fitensemble(data(trIdx,:),species(trIdx,:),'RUSBoost',20,t,'LearnRate',0.1);%See the documentation at:http://www.mathworks.com/help/stats/ensemble-methods.html#btgw1m1
                      yout = predict(B,data(teIdx,:));
                    
                      
                  case 'disc'
                      %Discriminant analysis. Assume each class has a Gaussian distribution with its own mean and covariance - It will give the model for each class
                      classstruct =  ClassificationDiscriminant.fit(data(trIdx,:),species(trIdx,:),'discrimType','diagQuadratic');
                      [yout score cost] = predict(classstruct,data(teIdx,:));
                            
                  case 'nn' %Nearest neighbot classifier
                            mdl = ClassificationKNN.fit(data(trIdx,:),species(trIdx,:),'NumNeighbors',3);
                            %L = loss(mdl,cancerdata(teIdx,:),cancerspecies(teIdx,:));
                            %disp(['K-NN lost: ' num2str(L)]);
                            yout = predict(mdl,data(teIdx,:)); %Predict the label of the output
      
                      
              end
              cp = classperf(species(teIdx,:),yout);
              cp_arr(i)=cp.CorrectRate;
              confusionmat(species(teIdx,:),yout);
              cfmat{i,1} = confusionmat(species(teIdx,:),yout);
              err(i) = sum(yout~=species(teIdx))/length(yout);
         end
        curcf = zeros(size(cfmat{1,1}));
        for cfmatidx=1:size(cfmat,1)
             curcf = curcf + cfmat{cfmatidx,1};
        end
        curcf = curcf/kval;
        cftotal = (cftotal*(roundidx-1)+curcf)/roundidx;
        %disp(['Round: ' num2str(roundidx)]);
        cftotalnew =cftotal./repmat(sum(cftotal,2),[1, size(cftotal,2)])*100;   
    end
end

function [feat,dataclass,samplename]=update_gland_texton_diagnosis_feature(filenames,labelpath,textonpath,mode,gradelist)
%Update the morphological features of all the glands
    ntextons = 50;
    feat = zeros(0,ntextons);
    dataclass = cell(0,1); %This is the class of the sample.
    ngrades = size(filenames,1);
    samplename = cell(0,1);
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
                %Load the texton histogram
                texton_idx_file_name = strcat(label_name,'_texton.mat'); %This is the label in large dataset
                load(strcat(textonpath,texton_idx_file_name),'new_text_map');
                glandidx = find(lblim==1);
                totalglandsize = length(glandidx);
                %Next, compute the rate for each band
                glandhist = hist(new_text_map(glandidx),[1:ntextons]);
               
                glandhist = glandhist/totalglandsize; %Normalization w.r.t gland size pixel sum.
                %Normalization w.r.t a single band
                %glandhist = glandhist/max(glandhist);
                dataclass(end+1,:)= gradelist(classidx);
                feat(end+1,:)=glandhist(:)';
                
                [meanNM,meanBPH,meanHGP,meanHGC,meanLGC,stdNM,stdBPH,stdHGP,stdHGC,stdLGC]=compute_feature_statistics(feat,dataclass);
                
                mean_arr = [meanNM;meanBPH;;meanLGC;meanHGP;meanHGC];
                std_arr = [stdNM;stdBPH;stdLGC;stdHGP;stdHGC];
                figure(2);
                barwitherr(std_arr(:,1:10)', mean_arr(:,1:10)');
                legend('NM','BPH','LGC','HGPIN','HGC');
                set(gca,'XTickLabel',{'T1','T2','T3','T4','T5','T6','T7','T8','T9','T10'})
                drawnow;
                figure(4);
                barwitherr(std_arr(:,11:20)', mean_arr(:,11:20)');
                legend('NM','BPH','LGC','HGPIN','HGC');
                set(gca,'XTickLabel',{'T11','T12','T13','T14','T15','T16','T17','T18','T19','T20'})
                figure(5);
                
                barwitherr(std_arr(:,21:30)', mean_arr(:,21:30)');
                legend('NM','BPH','LGC','HGPIN','HGC');
                set(gca,'XTickLabel',{'T21','T22','T23','T24','T25','T26','T27','T28','T29','T30'})
                drawnow;
                figure(6);
                barwitherr(std_arr(:,31:40)', mean_arr(:,31:40)');
                legend('NM','BPH','LGC','HGPIN','HGC');
                set(gca,'XTickLabel',{'T31','T32','T33','T34','T35','T36','T37','T38','T39','T40'})
                figure(7);
                barwitherr(std_arr(:,41:50)', mean_arr(:,41:50)');
                legend('NM','BPH','LGC','HGPIN','HGC');
                set(gca,'XTickLabel',{'T41','T42','T43','T44','T45','T46','T47','T48','T49','T50'})
               
                % drawnow;
                totalnumberofsampleprocessed = totalnumberofsampleprocessed + 1;
                disp(['Done...' num2str(totalnumberofsampleprocessed) ' samples']);
                
            
         
           end
    end
      
    
end

function [meanlt7,meaneq7,meangt7,stdlt7,stdeq7,stdgt7,scorelist]=compute_feature_statistics_cancer_class_on_score(featvect,dataclass)
    %Compute the mean and the standard deviation of all features for 5
    %different groups.
    %feat: a 2D array where each row is a data sample, each column is a
    %feature vector.
    meanfeat = mean(featvect,1);
    featvect = featvect./repmat(meanfeat,[size(featvect,1) 1]); %Normalize every feature so that every feature has mean value of 1
  
    nfeat = size(featvect,2);
    meanlt7 = zeros(1,nfeat);
    meaneq7 = zeros(1,nfeat);
    meangt7 = zeros(1,nfeat);
    stdlt7 = zeros(1,nfeat);
    stdeq7 = zeros(1,nfeat);
    stdgt7 = zeros(1,nfeat);
    scorelist = zeros(size(dataclass));
    idxlt7=[  find(strcmp(dataclass,'22')==1);....
                 find(strcmp(dataclass,'23')==1);...
                 find(strcmp(dataclass,'32')==1);...
                 find(strcmp(dataclass,'33')==1);...                 
              ];  
    if (~isempty(idxlt7))
        meanlt7=mean(featvect(idxlt7,:),1);
        stdlt7 =std(featvect(idxlt7,:),0,1);
        scorelist(idxlt7)=6;
    end
    
    idxeq7 = [find(strcmp(dataclass,'34')==1)
                 find(strcmp(dataclass,'43')==1)...
              ];  
    if (~isempty(idxeq7))
        meaneq7=mean(featvect(idxeq7,:),1);
        stdeq7 =std(featvect(idxeq7,:),0,1);
        scorelist(idxeq7)=7;
    
    end
    
    idxgt7 = [find(strcmp(dataclass,'44')==1);...
                 find(strcmp(dataclass,'45')==1);...
                 find(strcmp(dataclass,'45')==1);...
                 find(strcmp(dataclass,'55')==1);...
              ]; 
    if (~isempty(idxgt7))
        meangt7=mean(featvect(idxgt7,:),1);
        stdgt7 = std(featvect(idxgt7,:),0,1);
        scorelist(idxgt7)=8;
    
    end   
   
end

function [meanNM,meanHG,meanBP,meanCC,stdNM,stdHG,stdBP,stdCC,samplecond]=compute_feature_statistics_cancer_benign(featvect,dataclass)
    %Compute the mean and the standard deviation of all features for 5
    %different groups.
    %feat: a 2D array where each row is a data sample, each column is a
    %feature vector.
    meanfeat = mean(featvect,1);
    featvect = featvect./repmat(meanfeat,[size(featvect,1) 1]); %Normalize every feature so that every feature has mean value of 1
  
    nfeat = size(featvect,2);
    meanNM = zeros(1,nfeat);
    meanHG = zeros(1,nfeat);
    meanBP = zeros(1,nfeat);
    meanCC = zeros(1,nfeat); %For the cancer class
    stdNM = zeros(1,nfeat);
    stdHG = zeros(1,nfeat);
    stdBP = zeros(1,nfeat);
    stdCC = zeros(1,nfeat);
    samplecond = zeros(size(dataclass));%Condition
    idxNM=[  find(strcmp(dataclass,'NM')==1)                         
              ];  
    if (~isempty(idxNM))
        meanNM=mean(featvect(idxNM,:),1);
        stdNM =std(featvect(idxNM,:),0,1);
        samplecond(idxNM)=1;
    end
    
    idxHG = [find(strcmp(dataclass,'HG')==1)
              ];  
    if (~isempty(idxHG))
        meanHG=mean(featvect(idxHG,:),1);
        stdHG =std(featvect(idxHG,:),0,1);
        samplecond(idxHG)=2;
    
    end
    
    idxBP = [find(strcmp(dataclass,'BP')==1)]; 
    if (~isempty(idxBP))
        meanBP=mean(featvect(idxBP,:),1);
        stdBP = std(featvect(idxBP,:),0,1);
        samplecond(idxBP)=3;
    
    end   
   
    idxCC=[find(strcmp(dataclass,'23')==1);
           find(strcmp(dataclass,'33')==1);
           find(strcmp(dataclass,'32')==1);
           find(strcmp(dataclass,'34')==1);
           find(strcmp(dataclass,'43')==1);
           find(strcmp(dataclass,'44')==1);
           find(strcmp(dataclass,'45')==1);
           find(strcmp(dataclass,'54')==1);
              ];  
    if (~isempty(idxCC))
        meanCC=mean(featvect(idxCC,:),1);
        stdCC =std(featvect(idxCC,:),0,1);
        samplecond(idxCC)=4;
    end
end

function [hNM,hBPH,hHGP,hHGCancer,hLGCancer,x]=compute_histogram_of_features(featvect,dataclass)
    %Compute the histogram of the features vector for each class
    %Inputs:
    %   featvect: feature vector
    %   dataclass: class of each datasample
    %Outputs:
    %   hNM=>hLGCancer: histograms
    %   x: histogram bins
    minval = min(featvect);
    maxval = max(featvect);
    nbins = 50;
    x = linspace(minval,maxval,nbins);
    hNM = zeros(nbins,1);
    hBPH = zeros(nbins,1);
    hHGP = zeros(nbins,1);
    hHGCancer = zeros(nbins,1);
    hLGCancer = zeros(nbins,1);
    
    idxNM = find(strcmp(dataclass,'NM')==1);
    if (~isempty(idxNM))
        hNM=hist(featvect(idxNM),x);
        hNM=hNM/sum(hNM);
    end
    
    idxBP = find(strcmp(dataclass,'BP')==1);
    if (~isempty(idxBP))
        hBPH=hist(featvect(idxBP),x);
        hBPH = hBPH/sum(hBPH);
    end
    
    idxHGP = find(strcmp(dataclass,'HG')==1);
    if (~isempty(idxHGP))
        hHGP=hist(featvect(idxHGP),x);
        hHGP = hHGP/sum(hHGP);
    end
    
    idxLGCancer=[find(strcmp(dataclass,'23')==1);...
                 find(strcmp(dataclass,'32')==1);....
                 find(strcmp(dataclass,'33')==1)];
    if (~isempty(idxLGCancer))
        hLGCancer=hist(featvect(idxLGCancer),x);
        hLGCancer = hLGCancer/sum(hLGCancer);
    end
    
    idxHGCancer=[find(strcmp(dataclass,'35')==1);....
                 find(strcmp(dataclass,'44')==1);....
                 find(strcmp(dataclass,'45')==1);...
                 find(strcmp(dataclass,'53')==1);....
                 find(strcmp(dataclass,'54')==1);....
                 find(strcmp(dataclass,'55')==1);];
    if (~isempty(idxHGCancer))
        hHGCancer=hist(featvect(idxHGCancer),x);
        hHGCancer = hHGCancer/sum(hHGCancer);
    end
    
    
end
