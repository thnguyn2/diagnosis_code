function updateNuclei(filenames,curpath)
%This function computes the texton of single images and a final texton for
%all the images. The map of corresponding texton is also included
%odd and even symmetric filters
%Inputs:
%   filenames the name of all cores
%   curpath: path to the directory that saves response information information
%   datapath: path to the directory that saves the phase information.

%Outputs:
%   [None]

    addpath(curpath);
    
    ncols = 2048;
    nrows = 2048;
    ndirs = 20; %Number of pixel direction
    minrad = 8; %Radius of the voting region
    maxrad = 14;

    h = waitbar(0,'Overal percentage...');
    k = waitbar(0,'Hough voting...');
    
    %% Computing textons for each images..... 
    nfilestotal = 0;
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            nfilestotal =nfilestotal+1; 
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            fr_file_name = strcat(label_name,'_lm_fr_fine_scale.mat');%Name of the filter response
            hough_file_name = strcat(label_name,'_hough.mat');
            imname = strcat(label_name,'.mat');
            if (~exist(strcat(curpath,hough_file_name),'file'))
                disp(['Calculating nuclei for: ' label_name ' ...']);
                load(strcat(curpath,fr_file_name),'edir','emap','emag');
                edgemap = emap(:,:,1); %Edgemap
                dirmap = edir(:,:,1);  %Edge direction
                mask = im2bw(edgemap,0.5);
                s = regionprops(mask,'Area','PixelList');
                area_arrays = cat(1,s.Area);
                rejectidx = find(area_arrays<=3);%Get rid of the clutter
                newedgemap = edgemap;
                newdirmap = dirmap;
                for idx = 1:length(rejectidx)
                   curcoord=s(rejectidx(idx),1).PixelList;
                   coord1d = sub2ind([nrows ncols],curcoord(:,2),curcoord(:,1));
                   newedgemap(coord1d)=0;
                   newdirmap(coord1d)=0;
                end
                edgemap = newedgemap;
                dirmap = newdirmap;
                diredgemap = edgemap.*dirmap;
                diredgemap = (diredgemap-1)/ndirs*180; %(Angle here is in degree)
                idx = find(diredgemap<0); %Get rid of the vote of the background
                diredgemap(idx)=0;
                %Find all the pixel with non-zero index and process
                edgeidx = find(edgemap==1);
                radval = minrad:maxrad;
                radval = radval(:)';
                npixels = length(edgeidx);
                radval = repmat(radval,[npixels 2]);

                [y_coord,x_coord]=ind2sub([2048 2048], edgeidx); %Row and column indices of the edge pixels
                angleval = 180*(dirmap(edgeidx)-1)/ndirs;
                angleval = angleval(:);
                angleval = repmat(angleval,[1 size(radval,2)/2]); %Generate the first half of the angle
                angleval = [angleval (angleval+180)];
                %Now, compute the direction of the points that we should add

                x_offset = -sind(angleval).*radval;
                y_offset = cosd(angleval).*radval;
                nnbrs = size(x_offset,2);
                neigh_x =repmat(x_coord,[1 nnbrs]) + x_offset;
                neigh_y =repmat(y_coord,[1 nnbrs]) + y_offset;
                %Check for boundary conditions
                idx = find(neigh_x<1);
                neigh_x(idx)=1;
                idx = find(neigh_x>ncols);
                neigh_x(idx)=ncols;
                idx = find(neigh_y<1);
                neigh_y(idx)=1;
                idx = find(neigh_y>nrows);
                neigh_y(idx)=nrows;
                neigh_x = round(neigh_x);
                neigh_y = round(neigh_y);
                disp('Computing Hough map data...');
                houghmap = zeros(nrows,ncols);
                %Now, go for each neighbor and add up the count
                blob_std =5;
                hlog=-fspecial('log',[30 30],blob_std);
                
                for nbidx = 1:nnbrs
                    waitbar(nbidx/nnbrs,k,'Hough voting...')
                    votemask = zeros(nrows,ncols);
                    coord = [neigh_x(:,nbidx) neigh_y(:,nbidx)];
                    [coordval] = unique(coord,'rows');%Find the unique row elements
                    [temp1,temp2,idx] = unique([coordval;coord],'rows');%Find the unique row elements
                    npixelsassigned = size(coordval,1);
                    %Now, truncate the first part of the indexing so that only the
                    %overlapping's histogram is calculated
                    idx = idx(npixelsassigned+1:end);
                    nvotes = histc(idx,1:npixelsassigned);
                    updatedpixelidx = sub2ind([nrows ncols],coordval(:,2),coordval(:,1));
                    votemask(updatedpixelidx)=nvotes;
                    %Suppress false information at the edge
                    votemask(1,:)=0;
                    votemask(:,1)=0;
                    votemask(end,:)=0;
                    votemask(:,end)=0;
                    houghmap = houghmap + votemask;
                    figure(1);
                    imagesc(votemask);
                    colorbar
                    title('Vote map');
                    figure(2);
                    imagesc(houghmap);
                    colorbar;
                    title('Hough map');
                    drawnow;
                end
                figure(3)
                con_hough_map=imfilter(houghmap,hlog,'same');
                imagesc(con_hough_map);
                title('Concentrated Hough map');
                %Find the index of the center
                idx = find(con_hough_map>0.25*max(con_hough_map(:)));

                %Now, suppress the regions of blobs and replace it with the centrold
                mask = zeros(nrows,ncols);
                mask(idx)=1;
                %Before moving on, elminiate all regions with two few pixels
                mask = im2bw(mask,0.5);
                s = regionprops(mask,'Area','PixelList');
                area_arrays = cat(1,s.Area);
                rejectidx = find(area_arrays<=5);%Get rid of two small center
                for idx = 1:length(rejectidx)
                   curcoord=s(rejectidx(idx),1).PixelList;
                   coord1d = sub2ind([nrows ncols],curcoord(:,2),curcoord(:,1));
                   mask(coord1d)=0;
                end
                mask = im2bw(mask,0.5);
                s = regionprops(mask,'centroid');
                centroids = cat(1,s.Centroid);
                centroids = round(centroids);
                centroids1d = sub2ind([nrows ncols],centroids(:,2),centroids(:,1));
                %centroids1d = idx;
                edgemap(centroids1d)=2;
                figure(3)
                imagesc(edgemap);
                colorbar
                colormap gray
                nucleimap = zeros(nrows,ncols);
                nucleimap(centroids1d)=1;
                weightmap =zeros(nrows,ncols);
                weightmap(centroids1d)=con_hough_map(centroids1d);
                idx = find(nucleimap==0);
                minval = min(min((weightmap(centroids1d))));
                maxval = max(max((weightmap(centroids1d))));
                weightmap = 3*(weightmap-minval)/(maxval-minval)+1; %Normalize so that the min value is 2 and max is 3
                weightmap(idx)=0;
                save(strcat(curpath,hough_file_name),'nucleimap','houghmap','con_hough_map','weightmap','centroids1d','edgemap','-v7.3');
            else
                load(strcat(curpath,hough_file_name));
            end
            figure(4);
            imagesc(edgemap+weightmap);
            colorbar;
            title('Edge - Nuclei- Weight');
           
       end
    end
    close(k);
   
end