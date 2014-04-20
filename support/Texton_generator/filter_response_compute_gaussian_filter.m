function [f_res,conf_map,dir_map,scale_map]=filter_response_compute_gaussian_filter(im)
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
    dir_map=0;
    %Generate Leung Malik filter banks
    filters=makeGaussianFilters;
    
    nfilters=size(filters,3);
    nrows=size(im,1);
    ncols=size(im,2);
    %Compute the filter responses
    f_res=zeros(nrows,ncols,nfilters);
    disp('Filtering the images...');
    
    for filter_idx=1:nfilters
        f_res(:,:,filter_idx)=imfilter(im,filters(:,:,filter_idx),'same','circular');
    end
    
    nscales = 3;nang = 8;
    nrows = size(f_res,1);
    ncols = size(f_res,2);
    fe = zeros(nrows,ncols,nscales);
    fo = zeros(nrows,ncols,nscales);
     
    for scaleIdx=1:nscales
       curEvenBlock = f_res(:,:,(scaleIdx-1)*nang+1:scaleIdx*nang);
       curEvenBlock = curEvenBlock.^2;
       fe(:,:,scaleIdx)=sum(curEvenBlock,3);
    end
    
    %Compute the angle map by going through all possible angle and choose the one that gives maximum response. For each
    %angle, pick the maximum response by combining the odd and even filters
    angIm = zeros(nrows,ncols,nang);%Contains the maximum response for each direction....
    tempData = zeros(nrows,ncols,nscales);%Data contains information for each direction at multiple scales
    scaleData = zeros(nrows,ncols,nang);%This data contains the information on which scale gives the maximum response for each angle
    for angleIdx=1:nang
        for scaleIdx=1:nscales
            evenResIm = f_res(:,:,angleIdx+(scaleIdx-1)*nang);
            tempData(:,:,scaleIdx) = evenResIm.^2; %Compute the response at all scales
        end
        [angleIm(:,:,angleIdx),scaleData(:,:,angleIdx)] = max(tempData,[],3);
    end
    [tempMap,dir_map]=max(angleIm,[],3); %Produce a map for the angle
    for rowIdx=1:nrows
        for colIdx=1:ncols
            scale_map(rowIdx,colIdx)=scaleData(rowIdx,colIdx,dir_map(rowIdx,colIdx));%The Halo effect has to be treated
        end
    end
   conf_map = tempMap;
end
