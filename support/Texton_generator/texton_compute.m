function [texton,texton_map,filters,f_res,texton_diff]=texton_compute(im,k,f_type,disp_en,g_method,calc_en,fs_enable)
%Compute the texton and local image patches
%Written by: Tan H. Nguyen
%University of Illinois at Urbana-Champaign
%Data: Nov, 18th, 2012
%Version: 1.0
%Argument:
%   Inputs:
%       im: an input image
%       k: number of texton for grouping
%       disp_en: turn on/off displaying notifications
%       f_type: filter type ('lm' (Leung Malik), 's' (S), 'mr' (Maximum response))
%       g_method: grouping method ('kmean','kmed','greedy')
%       fs_enable: allow we do filtering on full-size images
%   Outputs:
%       textons: found texton response
%       local_im: local images corresponding to each texton
%       texton_map: a map of texton
%       filters: generated filters
%       f_res: filter responses
%       texton_diff: distance from the texton of each pixel to its central
%----------------------------------------------------------
    texton =0;
    texton_map=0;
    texton_diff=0;
    if (nargin<3)
        f_type = 'mr';
        g_method = 'kmean';
        disp_en=1;
    end
    if (nargin<7)
        fs_enable = 0;
    end
    if disp_en
        disp('Generating filters...');
    end
    if (fs_enable==0)
        switch f_type
            case {'lm'}
                 filters=makeLMfilters;
            case {'s'}
                filters=makeSfilters;
            case {'mr'}
                filters=makeRFSfilters;
        end
    else %Make filters with larger std.
         filters=makeLMfiltersFS;
    end
    nfilters=size(filters,3);
    nrows=size(im,1);
    ncols=size(im,2);
    %Compute the filter responses
    if (fs_enable==0)
            f_res=zeros(nrows,ncols,nfilters);
            if disp_en
                disp('Filtering the images...');
            end

            for filter_idx=1:nfilters
                %disp(['Working at scale... ' num2str(filter_idx)]);
                f_res(:,:,filter_idx)=imfilter(im,filters(:,:,filter_idx),'same','circular');
            end
    else
        %Compute the filtering for 10k x 10 k images goes here....
    end
    %If we have maximum response filter, just keep the maximum over all
    %direction
    if strcmp(f_type,'mr')
        new_res=zeros(nrows,ncols,8);
        for filter_idx=1:6
            new_res(:,:,filter_idx)=max(f_res(:,:,(filter_idx-1)*6+1:filter_idx*6),[],3);
        end
        new_res(:,:,7)=f_res(:,:,37);
        new_res(:,:,8)=f_res(:,:,38);
        f_res=new_res;
        clear new_res;
    end
    %Plot the filters if disp_en=1
    if (disp_en)
        sqrt_nfilt=ceil(sqrt(nfilters));
        figure;
        for filter_idx=1:nfilters
            subplot(sqrt_nfilt,sqrt_nfilt,filter_idx);
            imagesc(filters(:,:,filter_idx));
            
        end
        %colormap gray;
        title('Filter responses');
    end
    if (disp_en & calc_en)
        disp('Grouping and finding textons....');
        %Group the filter response to compute texton
        f_res_col_wise=zeros(ncols*nrows,size(f_res,3));
        for i=1:ncols
            for j=1:nrows
                f_res_col_wise((j-1)*ncols+i,:)=f_res(j,i,:);
            end
        end
        switch (g_method)
            case 'kmean'
                [texton_idx,texton]=kmeans_in_house(f_res_col_wise,k);%Do k-means on texton
            case 'kmed' %Not ready with large data yet. 
                %[texton_idx,texton]=kmedoids(f_res_col_wise,k);%Do k-means on texton
        end
        
        texton_map=zeros(nrows,ncols);
        texton_diff=zeros(nrows,ncols);
        %Create a texton_map
        for i=1:ncols
            for j=1:nrows
                texton_map(j,i)=texton_idx((j-1)*ncols+i);
                error_vect = squeeze(f_res(j,i,:))-texton(texton_map(j,i),:)';
                texton_diff(j,i)=norm(error_vect,'fro');
            end
        end
    end
    


end


function [local_patch]=compute_local_patch(filter_list,texton)
%   This function computes the local image patch given that it generates
%   the a texton (or a coef. vector) by filtering with multiple filters
    nfilters=size(filter_list,3); %Number of filters
    [nrows,ncols]=size(filter_list(:,:,1));
    F=zeros(nfilters,nrows*ncols);
    for i=1:nfilters
        F(i,:)=reshape(filter_list(:,:,i),[1 nrows*ncols]);
    end
    local_patch=F\texton;
    local_patch=reshape(local_patch,[nrows,ncols]);
end