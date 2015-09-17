function [classmap,stromamap,bgmap]=segmentation_refine(I,param)
%Refine the segmentation result using the watershed algorithm
    glandthresh=param.glandthresh; %Threshold for gland classification
    boundthreshforrejoin=param.boundthreshforrejoin;
    mindistancebetweenglands=param.mindistancebetweenglands;
    glandboundmin = param.glandboundmin;
    minlumenarea =param.minlumenarea;
    minglandsize = param.minglandsize;
    glandmap = (I<glandthresh);
    %First, get the image of the background
    lumenmap = (I==0);
    bgmapinv = imcomplement(lumenmap);
    bgmapinv = imclose(bgmapinv,strel('disk',40));
    bgmapinv = imfill(bgmapinv,'holes');
    bgmap = imcomplement(bgmapinv);
    bgidx = find(bgmap==1);
    lumenmap = lumenmap-bgmap;
    glandmap(bgidx)=0;%Clear the background
    glandmapinv = imcomplement(glandmap);
    glandmapinv = bwareaopen(glandmapinv, 2000);
    glandmap = imcomplement(glandmapinv);
    glandmap = imopen(glandmap,strel('disk',20));%Remove small connection
    glandmap = bwareaopen(glandmap,minglandsize);
    glandmapinv = imcomplement(glandmap);
    figure(1);imagesc(I);title('Original image');
    figure(2);imagesc(glandmap);title('Raw gland map estimation');
    glandmapdist = -bwdist(glandmapinv);
    regionmin = imextendedmin(glandmapdist,10);%Find all the regionmal min with metric distance to the boundary of at least 20
    %Suppress other local minima
    glandmapdist = imimposemin(glandmapdist,regionmin);
    figure(3);imshowpair(glandmap,regionmin,'blend');title('Gland map blended with regional minimum before watershed')
    figure(4);imagesc(glandmapdist);title('Distance map');colorbar;
    d=watershed(glandmapdist);
    boundaryidx = find(d==0);
    glandmap2 = glandmap;
    glandmap2(boundaryidx)=0;
    segmentcandidate = glandmap-glandmap2;
    minstromalengthforrejoin = 140;%Minimum distance for consider joining the parts
    segmentcandidate = imdilate(segmentcandidate,strel('disk',mindistancebetweenglands));%Make sure we have enough pixels to reliably evaluate the split
    edgepixidx = regionprops(im2bw(segmentcandidate),'PixelIdxList');
    ncandidate = size(edgepixidx,1);
    mapval = zeros(size(I));


    for idx=1:ncandidate
        curpixlist = edgepixidx(idx,1).PixelIdxList;
        glandmap2(curpixlist)=0;%Expand the gland boundary
        curmapval = mean(I(curpixlist));
        if ((curmapval<boundthreshforrejoin)&(length(curpixlist)>2*mindistancebetweenglands*minstromalengthforrejoin))
            glandmap2(curpixlist)=1; %Rejoin the gland
        end
        mapval(curpixlist)=length(curpixlist)/(2*mindistancebetweenglands);
    end
    cc = bwconncomp(glandmap2);l= labelmatrix(cc);rgbmap = label2rgb(l);
    figure(5);imagesc(rgbmap);title('Glandmap with boundary');
    figure(6);imagesc(mapval);colorbar; title('Mapval')
    segmentcandidate = glandmap-glandmap2;%Just add the new edge
    edgeidx = find(segmentcandidate==1);
    stromamap = (I>=0.5);
    stromamap = imclose(stromamap,strel('disk',10));
    classmap = 2*cast(stromamap,'single');
    glandmap2 = bwareaopen(glandmap2,minglandsize);
    glandidx=find(glandmap2==1);

    classmap(glandidx)=1;
    classmap(edgeidx)=2; %Turn the gland boundary into stroma

    %Now, go through each gland and get the lumen out
    glandboundingbox = regionprops(glandmap2,'BoundingBox');
    glandim = regionprops(glandmap2,'Image');

    nglands = size(glandboundingbox,1);

    for glandidx=1:nglands
         curbbox = glandboundingbox(glandidx,1).BoundingBox;
         c1 = round(curbbox(1));r1 = round(curbbox(2));
         c2 = c1+round(curbbox(3))-1;r2 = r1+round(curbbox(4))-1;
         curglandim = glandim(glandidx,1).Image;
         curlumenbox = lumenmap(r1:r2,c1:c2).*glandmap2(r1:r2,c1:c2);%Make sure we have all lumen inside the glands
         curlumenbox = bwareaopen(curlumenbox,minlumenarea);%Get rid too small lumen areas
         curlumenbox = imclose(curlumenbox,strel('disk',10));
         curlumenbox = imfill(curlumenbox,'holes');%Get rid of all debris
         curlumenbox = bwareaopen(curlumenbox,500).*curglandim;
         curlumenboxinv = imcomplement(curlumenbox);
         curlumenboxinv=bwareaopen(curlumenboxinv,minlumenarea);%Get rid of small areas like dots inside the lumen due to closing operation
         curlumenbox = imcomplement(curlumenboxinv);
         erodedim = imerode(curglandim,strel('disk',glandboundmin));
         erodedim(:,1:glandboundmin)=0;erodedim(:,end-glandboundmin:end)=0;%Clear all boundary as well
         erodedim(1:glandboundmin,:)=0;erodedim(end-glandboundmin:end,:)=0;
         curglandedge = curglandim - erodedim;
         newglandimage = (curglandim.*(~curlumenbox)|(curglandedge));%Preserve the gland boundary
         newglandimage = imclose(newglandimage,strel('disk',5));
         invnewglandimage = imcomplement(newglandimage);%Avoid bwareaopen on the background
         if (sum(newglandimage(:))>minlumenarea)%If the gland area is too small, do not run bwareaopen
            invnewglandimage = bwareaopen(invnewglandimage,minlumenarea);%Make sure we have no lumen area less than 20000 pixels
         end
         lumenpixidx = intersect(find(curglandim==1),find(invnewglandimage==1));
         [lumenpixrow,lumenpixcol]=ind2sub(size(curglandim),lumenpixidx); %Find local coordinates of lumen pixels
         lumenpixrow = r1+lumenpixrow-1;
         lumenpixcol = c1+lumenpixcol-1;
         lumenpixidx = sub2ind(size(I),lumenpixrow,lumenpixcol);
         classmap(lumenpixidx)=0;
         figure(9);imagesc(classmap);title('Final class map');
    end
    figure(8);imagesc(stromamap);title('Stromamap');
end