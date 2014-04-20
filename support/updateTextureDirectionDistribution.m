function updateTextureDirectionDistribution(filenames,labelpath,textonpath,histpath)
%This function computes the texton of single images and a final texton for
%all the images. The map of corresponding texton is also included
%odd and even symmetric filters
%Inputs:
%   filenames the name of all cores
%   labelpath,textonpath: paths to the label mask and texton indices
%   histpath: paths to the location for storing the histogram of directions
%Outputs:
    Fsym=makefilterfortexdir(1.0);
    Fnonsym= makefilterfortexdir(3.0);
    addpath(labelpath);
    addpath(textonpath)
    addpath('C:\Users\thnguyn2\Dropbox\current working code\Texton_generator');
    

    h = waitbar(0,'Calculating textons histogram...');
    
    %% Computing textons for each images.....    
    nrows = 2048;
    ncols = 2048;
    ndir = 8; 
    diameter = 40; %Radius of the area in which we compute the histogram of direction...
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
          
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            texture_dir_file_name = strcat(label_name,'_lm_fr.mat'); %For loading the texture direction
            texture_dir_hist_file_name =  strcat(label_name,'_text_dir_hist.mat');
            texton_idx_file_name = strcat(label_name,'_texton_hist_40.mat'); %For loading the label
            
            disp(['Processing image: ' label_name ' ...']);
            if (~exist(strcat(histpath,texture_dir_hist_file_name)))
                load(strcat(textonpath,texture_dir_file_name),'edir','emap');
                load(strcat(textonpath,texton_idx_file_name),'lblim');
                %Now, compute the histogram distribution of the texture
                %direction.
                nscales = 5; %
                hist = zeros(nrows,ncols,ndir);
                disp('Finding consistent over scales...');
                for idx = 1:ndir
                    finhist = ones(nrows,ncols);
                    for scaleidx=4:4
                        mask = zeros(nrows,ncols);
                        curim = edir(:,:,scaleidx);
                        pixidx = find(curim==idx); %Find satisfied pixels at current scale
                        mask(pixidx)=1;
                        finhist = finhist.*mask;
                    end
                    hist(:,:,idx)=finhist;
                end
                goodmask=sum(hist,3);
                newedir = edir(:,:,nscales).*goodmask;%Form an image of good texture direction from the largest scale
                figure(1);
                imagesc(newedir);
                title('Image of the best directions');
                colorbar
                                
                %% Working on this piece of code [1/23/2014]%%%
                disp('Calculating direction histogram (circular window)...');
                [rotvarhist,rotinvarhist]=computedirectionaldistribution(newedir, diameter);
                %%======================================================
                              
                %Find the entropy of the rotational variant histogram
                [maxval,maxidx]=max(rotvarhist,[],3);  
                rotvarhist = rotvarhist+1e-6;
                entrop=-log2(rotvarhist).*rotvarhist;
                entrop = sum(entrop,3);
                save(strcat(histpath,texture_dir_hist_file_name),'rotinvarhist',...
                    'maxval','maxidx','lblim','entrop');
               
            else
                load(strcat(histpath,texture_dir_hist_file_name),...
                    'maxval','maxidx','lblim','rotinvarhist','entrop');
                load(strcat(textonpath,texture_dir_file_name),'edir','emap');
            end
                     
             %This section of the code is for the paper...
             figure(1);
             imagesc(lblim);
             colorbar;
             truesize;
             title('Label image');
            
             figure(2);
             imagesc(rotinvarhist(:,:,1));
             colorbar;
             title('First bin of the histogram of direction...');
             figure(3);
             imagesc(entrop);
             axis off
             colormap jet;
             title('Entropy');
             colorbar;
             %End of the source code for the paper
            
                
       end
    end

    close(h);
end

function [rotvarhist,rotinvarhist,maxbinloc]=computedirectionaldistribution(edir, diameter)
%Compute the directional distribution at each pixels given the input image
%Inputs:
%   edir: a directional map for the texture direction at each pixel where
%   the background is set to 0
%   diameter: radius of the computing window
%Outputs:
%   rotvarhist,rotinvarhist: rotation variant and invariant histogram
%   distribution at each pixel...
    ndirs = 8;
    diameter = 40; 
    npixels = 2048^2;
    pixelidx=[1:npixels];
    pixelidx=pixelidx(:); %Line up current coordinates of the pixels
    x_offset = [-(diameter/2) (diameter/2) (diameter/2) -(diameter/2)];
    y_offset = [-(diameter/2) -(diameter/2) (diameter/2) (diameter/2)];
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
    rotvarhist = zeros(nrows,ncols,ndirs);
    rotinvarhist = zeros(nrows,ncols,ndirs);
    intIm = zeros(nrows,ncols,ndirs);
    disp('Calculating integral image....');
    for dirIdx=1:ndirs %Go through every texton
         idxIm = zeros(nrows,ncols);
         idx = find(edir==dirIdx);
         idxIm(idx)=1;
         idxIm = integralImage(idxIm);
         intIm(:,:,dirIdx) = idxIm(1:end-1,1:end-1);
    end
    
    %Compute the non-background count...
    intcountmap = integralImage(edir~=0);
    intcountmap = intcountmap(1:end-1,1:end-1);
    countmap = intcountmap(neigh_coord(:,3))+intcountmap(neigh_coord(:,1))-...
        intcountmap(neigh_coord(:,4))-intcountmap(neigh_coord(:,2));
     %Now, go through every texton again, for each texton, compute
     %the distribution pixel wise...
     for dirIdx = 1:ndirs
         curim = zeros(nrows,ncols);
         curintIm = intIm(:,:,dirIdx);
         curim(pixelidx) = curintIm(neigh_coord(:,3))+curintIm(neigh_coord(:,1))-...
         curintIm(neigh_coord(:,4))-curintIm(neigh_coord(:,2));
         %Normalize by how many pixels in the neighborhood
         curim(pixelidx)= curim(pixelidx)./countmap;
         rotvarhist(:,:,dirIdx)=curim;
     end
     %Rotate the histogram so that it has the order of reducing histogram
     %value
     rotvarhist2d = zeros(npixels,ndirs);
     rotinvarhist2d = zeros(npixels,ndirs);
     for dirIdx = 1:ndirs
         curhistmap = rotvarhist(:,:,dirIdx);
         rotvarhist2d(:,dirIdx)=curhistmap(:);
     end
     [maxval,maxloc]=max(rotvarhist2d,[],2);
     %Perform rotational shifting...
     rowidx = [1:npixels];
     rowidx=rowidx(:);
     for orderIdx=1:ndirs
         curcolidx = mod(maxloc(:)+orderIdx-1,ndirs);%This is the current order of the orderIdx-th max elements
         curcolidx(find(curcolidx==0))=ndirs; %Make sure that the column index range is from 1 to ndirs
         rotinvarhist2d(:,orderIdx)=rotvarhist2d((curcolidx-1)*npixels+rowidx);
         %Convert back into 2d map
         rotinvarhist(:,:,orderIdx)=reshape(rotinvarhist2d(:,orderIdx),[nrows ncols]);
     end     
     maxbinloc = reshape(maxloc,[nrows ncols]);
end

function F=makefilterfortexdir(axisscale)
% Returns the LML filter bank of size 49x49x48 in F. To convolve an
% image I with the filter bank you can either use the matlab function
% conv2, i.e. responses(:,:,i)=conv2(I,F(:,:,i),'valid'), or use the
% Fourier transform.

  SUP=60;                 % Support of the largest filter (must be odd)
  SCALEX=30;  % Sigma_{x} for the oriented filters
  NORIENT=8;              % Number of orientations

  NF = NORIENT;
  F=zeros(SUP,SUP,NF);
  hsup=(SUP-1)/2;
  [x,y]=meshgrid([-hsup:hsup],[hsup:-1:-hsup]);
  orgpts=[x(:) y(:)]';

  count=1;
  for orient=0:NORIENT-1,
      angle=pi*orient/NORIENT;  % Not 2pi as filters have symmetry
      c=cos(angle);s=sin(angle);
      rotpts=[c -s;s c]*orgpts;
      F(:,:,count)=makefilter(SCALEX,0,0,rotpts,SUP,axisscale);    %Take the first derivative in x and the second derivative in y
      count=count+1;
  end;
end

function f=makefilter(scale,phasex,phasey,pts,sup,axisscale)
%Phase x, phase y is the order of Gaussian function
  gx=gauss1d(scale,0,pts(1,:),phasex);
  gy=gauss1d(scale/axisscale,0,pts(2,:),phasey);
  f=normalise(reshape(gx.*gy,sup,sup));
end

function g=gauss1d(sigma,mean,x,ord)
% Function to compute gaussian derivatives of order 0 <= ord < 3
% evaluated at x.
% ord is the order of gaussian derivative

  x=x-mean;num=x.*x;
  variance=sigma^2;
  denom=2*variance;  
  g=exp(-num/denom)/(pi*denom)^0.5;
  switch ord,
      case 0, g=g;
    case 1, g=-g.*(x/variance);                 %First order derivative of Gaussian
    case 2, g=g.*((num-variance)/(variance^2)); %Second order derivative of gaussian
  end;
end

function f=normalise(f)
    f=f/sum(abs(f(:)))
end