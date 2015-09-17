%                     cribriformarr=zeros(nregions,1);
%                     %Compute the compactness of different glands...
%                     glandcompactness = zeros(nregions,1);
%                     minlumenareaforcribriform = 1500;
%                     cribriformeccentthresh = 0.85;
%                     %glandmapwithconvexgland=zeros(sizelblim); %This is the map of all glandmap.
%                     curF=zeros(nregions,P);%Add the new feature to the feature map
%                     for regionidx=1:nregions
%                     %for regionidx=8:8
%                         
%                         %Compute the cribriform score by computing the
%                         %number of lumen area inside each gland. The
%                         %surrounded lumen is contrainst to have two main
%                         %axis of ratio from 0.8 to 1.2
%                         curmask = zeros(size(lblim));
%                         curmask(pixellistarr(regionidx).PixelIdxList)=1;
%                         curmask = curmask(r1:r2,c1:c2);
%                             
%                           
%                         tempim = lblim(r1:r2,c1:c2).*curmask;
%                         tempglandmap = im2double((tempim==1));
%                         %Next, close the image, get its boundary and
%                         %close the gland
%                         closingtempglandmap = imclose(tempglandmap,strel('disk',20));
%                         %calculating the perimeter of the closed glandmap
%                         filledclosedtempmap=imfill(closingtempglandmap,'hole');
%                         filledglandperi=regionprops(filledclosedtempmap,'Perimeter');
%                         boundaryperiratioarr(regionidx,1)=filledglandperi(1).Perimeter/convexhullperi(1).Perimeter;
%                         
%                         glandsedge = closingtempglandmap - imerode(closingtempglandmap,strel('disk',3)); 
%                         %Closing the image so that it brings the gap
%                         %between small area
%                         edgeidx = find(glandsedge==1);
%                         tempglandmap(edgeidx)=1;%Create a pseudo edge in the region. This is good for regions with very thin gap
%                         tempglandmapnohole = imfill(tempglandmap,'holes');
%                             [lumenrowidx,lumencolidx] = find(holemap==1); %Display all the lumen area
%                             lumenrowidx=r1+lumenrowidx-1;
%                             lumencolidx=c1+lumencolidx-1;
%                             lblim((lumencolidx-1)*nrows+lumenrowidx)=3;
%                             
%                             lumeneccent = regionprops(holemap,'Eccentricity'); %Get the eccentricity of the lumen area
%                            
%                             %Compute the cribriform score by counting the
%                             %number of round area inside each gland.
%                             if (isempty(lumeneccent))
%                                 cribriformarr(regionidx,1)=0;
%                             else
%                                 nlumens = size(lumeneccent,1);%Count the number of lumens
%                                 lumeneccentarr = zeros(nlumens,1);
%                                 for lumenidx = 1:nlumens
%                                     lumeneccentarr(lumenidx,1)=lumeneccent(lumenidx).Eccentricity;
%                                 end
%                                 [goodcribriformlumenidx]=find(lumeneccentarr<cribriformeccentthresh);
%                                 if (~isempty(goodcribriformlumenidx))
%                                     cribriformarr(regionidx,1) = abs(length(goodcribriformlumenidx)); %1 is for eliminating the only lumen area
%                                 else
%                                     cribriformarr(regionidx,1) = 0;
%                                 end
%                             end
%                         end
%                     
%                     %Calcuate the lumen parameters...
%                     lumenarea = zeros(nregions,1); 
%                     nonsamplemap = (lblim~=0);
%                     nonsamplemap = imclose(nonsamplemap,strel('disk',30)); %Avoid the phenomenon when the lumen is cut through
%                     nonsamplemap = imfill(nonsamplemap,'holes');
% 
%                     lumenmap = nonsamplemap.*(lblim==0); %Get lumens within convex hull
%                     %Clear out small lumen area at the gland's boundary
%                     lumenmap = imopen(lumenmap,strel('disk',10));
%                     lumenmap = imfill(lumenmap,'holes');
%                     %Remove areas of lumen with too few pixels
%                     lumenmap = bwareaopen(lumenmap,2000,8);
% 
%                     %Compute the lumen area for each gland...
%                     pixidxperregion = regionprops(glandmap,'PixelIdxList');
%                     glandeccentarr = zeros(nregions,1);
%                     figure(3);
%                     imagesc(lblim);
%                     for glandregionidx = 1:nregions
%                          lumenarea(glandregionidx,1) = sum(lumenmap(pixidxperregion(glandregionidx,1).PixelIdxList));
%                           
%                     end
%                     
%                     curF(:,1) = glandarea/mean(glandarea);%First, we normalize every feature w.r.t to the mean
%                   
%                     drawnow;
%                     lumentoglandratioarr = lumenarea./glandarea;
%                     %% First feature: Standard deviation of the gland size/its mean
%                     glandstdovermeanarea = std(glandarea)/mean(glandarea);
%                     
%                     %% Second feature: Gland's eccentricity                    
%                     [mugldecc,stdglecc,skewglecc,kurtval,maxglecc]=get_pdf_moment(glandeccentarr,false,glandarea);
%                     %Compute the mean scallop ratio
%                     %scallopratiomean = mean(scallopratioarr.*glandarea)/mean(glandarea);
%                                         
%                     [muperi,stdperi,skewperi,kurtperi,maxperi]=get_pdf_moment(boundaryperiratioarr,false,glandarea);
%                     [muscal,stdscal,skewscal,kurtscal]=get_pdf_moment(scallopratioarr,true,glandarea);
%                     [musi,stdsi,skewsi,kurtsi,maxsi]=get_pdf_moment(stromainvasionarray,true,glandarea);
%     %                periratiomean = mean(boundaryperiratioarr.*glandarea)/mean(glandarea);
%                     [mucri,stdcri,skewcri,kurtcri,maxcri]=get_pdf_moment(cribriformarr,true,glandarea);
%                      

%                     feat(end+1,:)= [glandstdovermeanarea mugldecc skewglecc maxglecc muscal skewscal muperi skewperi maxperi musi skewsi maxsi mucri];
%                     indivglandfeat{end+1,1}=curF;
%                     
%                     dataclass(end+1,:)= gradelist(classidx);
% 
%                     %Plot the feature distribution with a colorbar
%                     drawnow;
%                     figure(1);
%                     barwitherr(std_arr(:,13)', mean_arr(:,13)');
%                     title('Perimeter ratio score')
%                     disp(['Current feature: Mean ' num2str(muperi) ', Min: ' num2str(min(boundaryperiratioarr)) ', Max:' num2str(max(boundaryperiratioarr))]);
%                    
%                     
%                     totalnumberofsampleprocessed = totalnumberofsampleprocessed + 1;
%                    % disp(['Done...' num2str(totalnumberofsampleprocessed) ' samples']);
