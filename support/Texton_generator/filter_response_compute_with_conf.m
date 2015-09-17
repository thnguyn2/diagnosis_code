function [f_res,mage,edgebordere,edgebordero,edgeborder,dir_map]=filter_response_compute_with_conf(im)
%Compute the confidence map, normal filter response and the gradient
%direction map
%Written by: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
%Data: Nov, 18th, 2012
%Version: 1.0
%Argument:
%   Inputs:
%       im: an input image
%   Outputs:
%       f_res: filter response
%       conf_map: confidence map
%       dir_map: map of the gradient directon
%       scale_map: map of the scale which gives maximum resposne
%----------------------------------------------------------
    thresh = 0.2;
    %Generate Leung Malik filter banks
    filters=makeLMfilters;
      
    nfilters=size(filters,3);
    nrows=size(im,1);
    ncols=size(im,2);
    %Compute the filter responses
    f_res=zeros(nrows,ncols,nfilters);
    disp('Filtering the images...');
    
    for filter_idx=1:nfilters
        f_res(:,:,filter_idx)=imfilter(im,filters(:,:,filter_idx),'same','circular');
    end
    
    nscales = 5;nang = 8; ninvscales = 6;
    nrows = size(f_res,1);
    ncols = size(f_res,2);
    
     
    %Go through all scales. For each scale, go with the compute the maximum
    %response and choose the best angle
    edgeborder = zeros(nrows,ncols,nscales);
    edgebordero = zeros(nrows,ncols,nscales);
    edgebordere = zeros(nrows,ncols,nscales);
    mag=zeros(nrows,ncols,nscales);
    mage =zeros(nrows,ncols,nscales);
    mago =zeros(nrows,ncols,nscales);
    curEvenData=zeros(nrows,ncols,nang);
    curOddData=zeros(nrows,ncols,nang);
    dir_map=zeros(nrows,ncols,nscales);
   

    basesize=nrows*ncols;
    id=[1:basesize]';
    for scaleIdx=1:nscales
        for angleIdx=1:nang
             curOddData(:,:,angleIdx) = f_res(:,:,angleIdx+(scaleIdx-1)*nang);
            curEvenData(:,:,angleIdx) = f_res(:,:,angleIdx+(scaleIdx-1)*nang+nang*nscales);
        end
        curScaleData=sqrt(curOddData.^2+curEvenData.^2);
        [temp,maxIdx]=max(curScaleData,[],3); %Search for the gradient angle by choosing the max direction over all posible angles
        % of the directional gradient
        %Extract the even part
        mage(:,:,scaleIdx) = reshape(curEvenData(id+(maxIdx(:)-1)*basesize),[nrows,ncols]);%What does this one do?
        mag(:,:,scaleIdx)=sqrt(sum(curOddData.^2,3)+sum(curEvenData.^2,3)); %This line will help avoiding the Winner takes all clase
        mago(:,:,scaleIdx) = reshape(curOddData(id+(maxIdx(:)-1)*basesize),[nrows,ncols]);
        curmag = mag(:,:,scaleIdx);
        curmage = mage(:,:,scaleIdx);
        curmago = mago(:,:,scaleIdx);
        %Find the positive location
        curmage = (curmage>0)-(curmage<0); %Use a historisis thresholding here...
        maxmagoval =max(abs(curmago(:)));
        %curmago = (curmago>0.1*maxmagoval)-(curmago<-0.1*maxmagoval);
        curmago = (curmago>0)-(curmago<0);

        %Find the zero-crossing
        hidx=[(curmage(2:nrows,:)~=curmage(1:nrows-1,:));zeros(1,ncols)];
        vidx=[(curmage(:,2:ncols)~=curmage(:,1:ncols-1)),zeros(nrows,1)];
        
        %Find the zero-crossing for the odd image
        hoidx=[(curmago(2:nrows,:)~=curmago(1:nrows-1,:));zeros(1,ncols)];
        voidx=[(curmago(:,2:ncols)~=curmago(:,1:ncols-1)),zeros(nrows,1)];
        
        edgebordere(:,:,scaleIdx)=(hidx|vidx).*(curmag>thresh*max(curmag(:)));
        edgebordero(:,:,scaleIdx)=(hoidx|voidx).*(curmag>thresh*max(curmag(:)));
        edgeborder(:,:,scaleIdx) = edgebordero(:,:,scaleIdx)+ 2*edgebordere(:,:,scaleIdx);
        dir_map(:,:,scaleIdx)  = (maxIdx-1)*180/nang;
    end
    
    %Go through all scales. For each scale, go with the compute the maximum
    %response and choose the best angle
    sedgeborder = zeros(nrows,ncols,ninvscales);
    smag=zeros(nrows,ncols,ninvscales);
    sfres = zeros(nrows,ncols,ninvscales);
    for scaleIdx=1:ninvscales
        sfres(:,:,scaleIdx) = f_res(:,:,2*nang*nscales+2*scaleIdx);
        curres=sfres(:,:,scaleIdx);
        curres = (curres>0)-(curres<0);
        %Find the zero-crossing for the odd image
        shoidx=[(curres(2:nrows,:)~=curres(1:nrows-1,:));zeros(1,ncols)];
        svoidx=[(curres(:,2:ncols)~=curres(:,1:ncols-1)),zeros(nrows,1)];
        sedgeborder(:,:,scaleIdx) = (shoidx|svoidx).*(curmag>thresh*max(curmag(:)));;
    end
    %Now, compute the filter response based on norm
    
end
