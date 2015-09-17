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
    ntextons = 50;
    retrain = 0;
    h = waitbar(0,'Calculating textons...');
    radius = 60;%In the 3072 x 3072 pixels image
   
    nrows = 3072;
    ncols = 3072;
 
    npixels = nrows*ncols;
    
    pixelidx=[1:npixels];
    pixelidx=pixelidx(:); %Line up current coordinates of the pixels
    x_offset = [-radius radius radius -radius];
    y_offset = [-radius -radius radius radius];
    [y_coord,x_coord]=ind2sub([nrows ncols], pixelidx); %Row and column indices of current pixels
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
    %ngrades = size(filenames,1);
    nsamples = size(filenames,1);
    for sampleIdx=1:nsamples
    %for classidx=1:ngrades %Go through different classes
    %   nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
    %   for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
          
           % cur_file_name = filenames{classidx,1}{sampleIdx,1};
            cur_file_name = filenames{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            texton_idx_file_name = strcat(label_name,'_texton.mat'); %This is the label in large dataset
            texton_hist_file_name =  strcat(label_name,'_texton_hist_',num2str(radius),'.mat');
            
            disp(['Processing image: ' label_name ' ...']);
            %lbl_file_name = strcat(label_name,'_resized.mat');
            computingPlatform = 0; 
            if ((~exist(strcat(textonpath,texton_hist_file_name)))||(retrain==1))
                load(strcat(textonpath,texton_idx_file_name),'new_text_map');
%                load(strcat(labelpath,lbl_file_name));
                %Downsample the texton index image and the lbl image to
                %save some space
                new_text_map = imresize(new_text_map,[nrows ncols],'nearest');
%               lblim = imresize(lblim,[nrows ncols],'nearest');
                %Convert into the pixel-wise distribution
                histim = zeros(nrows,ncols,ntextons);
                intIm = zeros(nrows,ncols,ntextons);
                intIm = cast(intIm,'single');
                disp('Calculating integral image....');
                hist_filter = fspecial('gaussian',[round(2*radius) round(2*radius)],radius/5);
                tic;
                
                for textonIdx=1:ntextons %Go through every texton
                   idxIm = zeros(nrows,ncols);
                   idx = find(new_text_map==textonIdx);
                   idxIm(idx)=1;
                   if (computingPlatform==0)
                        idxIm = imfilter(idxIm,hist_filter,'same');
                   else
                   end
                   idxIm = integralImage(idxIm);
                   intIm(:,:,textonIdx) = cast(idxIm(1:end-1,1:end-1),'single');
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
                histim = cast(histim,'single');
                if (exist('lblim'))
                    save(strcat(textonpath,texton_hist_file_name),'histim','lblim','-v6');
                else
                    save(strcat(textonpath,texton_hist_file_name),'histim','-v6');
                
                end
%                 figure(1);
%                 imagesc(lblim);
%                 title('Label map');
%                 drawnow;

                for textonIdx=1:20:20
                     figure(textonIdx+1);
                     imagesc(histim(:,:,textonIdx));
                     colorbar;
                end
                timeElap = toc;
                disp(['Time elapse: ' num2str(timeElap)]);
               
            else
               % load(strcat(textonpath,texton_hist_file_name),'histim','lblim');
                
            end
          
           
            
            
     %  end
    end

    close(h);
end