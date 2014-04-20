function updateTextonDistribution(filenames,labelpath,textonpath)
%This function computes the texton of single images and a final texton for
%all the images. The map of corresponding texton is also included
%odd and even symmetric filters
%Inputs:
%   filenames the name of all cores
%   labelpath,textonpath: paths to the label mask and texton indices
%Outputs:
    
    addpath(labelpath);
    addpath(textonpath)
    addpath('C:\Users\thnguyn2\Dropbox\current working code\Texton_generator');
    ntextons = 51;

    h = waitbar(0,'Calculating textons...');
    %radius = 20;%The window size is 32x4 in the real image, (about 128/14 approx 9 microns^2 window)
   
    radius = 40;%The window size is 32x4 in the real image, (about 128/14 approx 9 microns^2 window)
    
    npixels = 2048^2;
    pixelidx=[1:npixels];
    pixelidx=pixelidx(:); %Line up current coordinates of the pixels
    x_offset = [-(radius/2) (radius/2) (radius/2) -(radius/2)];
    y_offset = [-(radius/2) -(radius/2) (radius/2) (radius/2)];
    nrows = 2048;
    ncols = 2048;
    [y_coord,x_coord]=ind2sub([2048 2048], pixelidx); %Row and column indices of current pixels
    y_coord = y_coord(:);
    x_coord = x_coord(:);
    diff_coord = x_offset*nrows+y_offset; %Difference to the center pixel
    %Calculate the coordinates of 4 corners of all pixels
    neigh_x = repmat(x_coord,[1 4])+ones(npixels,1)*x_offset;
    neigh_y = repmat(y_coord,[1 4])+ones(npixels,1)*y_offset;
    %Check for boundary conditions
    idx = find(neigh_x<1);
    neigh_x(idx)=1;
    idx = find(neigh_x>ncols);
    neigh_x(idx)=ncols;
    idx = find(neigh_y<1);
    neigh_y(idx)=1;
    idx = find(neigh_y>nrows);
    neigh_y(idx)=nrows;
    %Convert into linear indexing....
    neigh_coord=sub2ind([nrows ncols],neigh_y,neigh_x);
    %Compute the number of neighbors for each pixels
    nneighbors = (neigh_x(:,2)-neigh_x(:,1)+1).*(neigh_y(:,3)-neigh_y(:,2)+1);
    %% Computing textons for each images..... 
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
          
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            texton_idx_file_name = strcat(label_name,'_texton.mat'); %This is the label in large dataset
            texton_hist_file_name =  strcat(label_name,'_texton_hist_',num2str(radius),'.mat');
            
            disp(['Processing image: ' label_name ' ...']);
            lbl_file_name = strcat(label_name,'.mat');

            if (~exist(strcat(textonpath,texton_hist_file_name)))
                load(strcat(textonpath,texton_idx_file_name),'new_text_map','ovr_posterior');
                load(strcat(labelpath,lbl_file_name));
                %Downsample the images to match the size
                lblimnew = imresize(lblim,[2048 2048],'nearest');
               
                lblim = lblimnew;%Resize the image to match the same size of the label
 
                %Convert into the pixel-wise distribution
                histim = zeros(nrows,ncols,ntextons);
                intIm = zeros(nrows,ncols,ntextons);
                disp('Calculating integral image....');
                for textonIdx=1:ntextons %Go through every texton
                    idxIm = zeros(nrows,ncols);
                    idx = find(new_text_map==textonIdx);
                    idxIm(idx)=1;
                    idxIm = integralImage(idxIm);
                    intIm(:,:,textonIdx) = idxIm(1:end-1,1:end-1);
                end

                disp('Calculating histgram of texton....');
                %Now, go through every texton again, for each texton, compute
                %the distribution pixel wise...
                for textonIdx = 1:ntextons
                    curim = zeros(nrows,ncols);
                    curintIm = intIm(:,:,textonIdx);
                    curim(pixelidx) = curintIm(neigh_coord(:,3))+curintIm(neigh_coord(:,1))-...
                        curintIm(neigh_coord(:,4))-curintIm(neigh_coord(:,2));
                    %Normalize by how many pixels in the neighborhood
                    curim(pixelidx)= curim(pixelidx)./nneighbors;
                    histim(:,:,textonIdx)=curim;
                end
                save(strcat(textonpath,texton_hist_file_name),'histim','lblim','-v7.3');
            else
                load(strcat(textonpath,texton_hist_file_name),'histim','lblim');
                
            end
            load(strcat(textonpath,texton_idx_file_name),'ovr_posterior');
            figure(1);
            imagesc(lblim);
            title('Label map');
            drawnow;
            
            for textonIdx=1:17:ntextons
                 figure(textonIdx+1);
                 imagesc(histim(:,:,textonIdx));
                 colorbar;
            end
            figure(ntextons+1);
            imagesc(ovr_posterior);
            colorbar;
            title('Posterior');          
            
            
            
       end
    end

    close(h);
end