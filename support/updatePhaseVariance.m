function updatePhaseVariance(filenames,labelpath,phasevarpath)
%This function computes the texton of single images and a final texton for
%all the images. The map of corresponding texton is also included
%odd and even symmetric filters
%Inputs:
%   filenames the name of all cores
%   labelpath,phasevarpath: paths to the label mask and the location for phase variance
%Outputs:
    
    addpath(labelpath);
    addpath(phasevarpath)
  
    h = waitbar(0,'Calculating phase variance...');
    radius = 64;%The window size is 32x4 in the real image, (about 128/14 approx 9 microns^2 window)
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

   %% Computing textons for each images..... 
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            nneighbors = (neigh_x(:,2)-neigh_x(:,1)+1).*(neigh_y(:,3)-neigh_y(:,2)+1);
            waitbar(sampleIdx/nsamples,h,'Progress...')
          
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            label_file_name = strcat(label_name,'.mat');
            phase_var_file_name = strcat(label_name,'_phase_var.mat'); %This is the label in large dataset
            disp(['Processing image: ' label_name ' ...']);
        
            if (~exist(strcat(phasevarpath,phase_var_file_name)))
                org_im = imread(cur_file_name); %Load the original image
                load(strcat(labelpath,label_file_name));
                
                %Resize the image
                newim = imresize(org_im,[nrows,ncols]);
                lblimnew = imresize(lblim,[nrows, ncols],'nearest');
                newim = cast(newim,'single');
                %First, compute the mean phase phase value value in a window
                disp('Calculating mean value....Halo included');
                intIm = integralImage(newim);
                intIm = intIm(1:end-1,1:end-1);
                meanim = zeros(nrows,ncols);
                meanim(pixelidx) = intIm(neigh_coord(:,3))+intIm(neigh_coord(:,1))-...
                        intIm(neigh_coord(:,4))-intIm(neigh_coord(:,2));
                meanim(pixelidx)=meanim(pixelidx)./nneighbors;
                meanim = cast(meanim,'single');
                
                disp('Calculating phase variance....Halo included');
                diff = newim - meanim;
                diff = diff.^2;
                intIm = integralImage(diff);
                intIm = intIm(1:end-1,1:end-1);
                varim = zeros(nrows,ncols);
                varim(pixelidx) = intIm(neigh_coord(:,3))+intIm(neigh_coord(:,1))-...
                        intIm(neigh_coord(:,4))-intIm(neigh_coord(:,2));
                varim(pixelidx)=varim(pixelidx)./nneighbors;
                stdim = sqrt(varim);
                
                %Now, compute the mean and phase variance in the region
                %without including Halo
                nonhaloidx = find(lblimnew~=0);
                zeroidx = find(lblimnew==0);
                %First, zero-out the contribution of the pixel with Halo or
                %lumen.
                newim(zeroidx)=0;
                %Next, compute the mean phase
                calcMask = zeros(nrows,ncols);
                calcMask(nonhaloidx)=1;
                nneighbors = zeros(nrows,ncols);
                disp('Compute number of neighbors...Halo rejected');
                intIm = integralImage(calcMask);
                intIm = intIm(1:end-1,1:end-1);
                nneighbors(pixelidx) = intIm(neigh_coord(:,3))+intIm(neigh_coord(:,1))-...
                        intIm(neigh_coord(:,4))-intIm(neigh_coord(:,2)); %Number of good neighbors
                
                disp('Calculating mean value....Halo rejected');
                intIm = integralImage(newim);
                intIm = intIm(1:end-1,1:end-1);
                meanimnoHalo = zeros(nrows,ncols);
                meanimnoHalo(pixelidx) = intIm(neigh_coord(:,3))+intIm(neigh_coord(:,1))-...
                        intIm(neigh_coord(:,4))-intIm(neigh_coord(:,2));
                meanimnoHalo=meanimnoHalo./(nneighbors+1e-3);
                meanimnoHalo(zeroidx)=0;
                meanimnoHalo = cast(meanimnoHalo,'single');                
                 
                meanimnoHalozeroed =  meanimnoHalo;
                meanimnoHalozeroed(zeroidx)=0;
                newim = cast(newim,'single');
                disp('Calculating phase variance....Halo rejected');
                diff = newim -meanimnoHalozeroed;
                diff = diff.^2;
                intIm = integralImage(diff);
                intIm = intIm(1:end-1,1:end-1);
                varimnoHalo = zeros(nrows,ncols);
                varimnoHalo(pixelidx) = intIm(neigh_coord(:,3))+intIm(neigh_coord(:,1))-...
                        intIm(neigh_coord(:,4))-intIm(neigh_coord(:,2));
                varimnoHalo=varimnoHalo./nneighbors;
                stdimnoHalo = sqrt(varimnoHalo);
                stdimnoHalo(zeroidx)=0;              
                
                save(strcat(phasevarpath,phase_var_file_name),'meanim','meanimnoHalo',...
                    'stdim','stdimnoHalo','lblim','lblimnew');
            else
                load(strcat(phasevarpath,phase_var_file_name),'meanim','meanimnoHalo','stdim',...
                    'stdimnoHalo','lblimnew');
              
            end
            pause;
            figure(1);
            meanim = meanim/65536*4.0-0.5;
            imagesc(meanim);
            colorbar
            title('Mean Halo');
            figure(2)
            stdim = stdim/65536*4-0.5;
            imagesc(stdim);
            colorbar
            title('Std Halo');
            figure(3)
            imagesc(lblimnew);
           
            drawnow;
            
            figure(4);
            meanimnoHalo = meanimnoHalo/65536*4.0-0.5;
            meanimnoHalo = medfilt2(meanimnoHalo,[3 3]);
            imagesc(meanimnoHalo);
            colorbar
            title('Mean no Halo');
            figure(5)
            stdimnoHalo = stdimnoHalo/65536*4-0.5;
            stdimnoHalo = medfilt2(stdimnoHalo,[3 3]);
            imagesc(stdimnoHalo);
            colorbar
            title('Std no Halo');
            
            figure(6)
            imagesc(lblimnew);
            
           % pause;
            
       end
    end

    close(h);
end