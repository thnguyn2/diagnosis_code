function updateLs(filenames,curpath)
%This function computes ls for different images
%Inputs:
%   filenames, glandnames: the name of all cores
%   path: the place to save phase distributions
%Outputs:
%   dist: a 4x1 cell array containing distributions corresponding to 4 different
%   classes. Each cell is a 3x500 matrix where each row is a histogram
%   distribution of the phase value...
    
    %part1: compute the distribution for each image and save it
    addpath(curpath);
    h = waitbar(0,'Calculating scattering mean free path...');
    hy = fspecial('sobel');
    hx = hy';
    hs = fspecial('gaussian',[3 3],0.5); %Filter the image first before computing the gradient
    blk_size = [20 20]; %Size that we cmpute the variance
    newsize = 2048;     %Size of the new image
    
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       lsim = zeros(newsize,newsize);
       gim = zeros(newsize,newsize);
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            ls_file_name = strcat(label_name,'_ls.mat');
            if (~exist(ls_file_name,'file'))
                disp(['Calculating ls: ' label_name ' ...']);
                org_im = imread(cur_file_name); %Load the original data
                %Compute the gradient of the orginal image
                org_im = im2double(org_im);
                org_im = phase_im_convert(org_im); %Convert the image into phase radian values
                
                %Now, resize the image so that its dimension is a multiple
                %of that of the downsample image
                nrows = size(org_im,1);
                ncols = size(org_im,2);
                newnrows = ceil(nrows/newsize)*newsize;
                newncols = ceil(nrows/newsize)*newsize;
                org_im = imresize(org_im,[newnrows,newncols]);
                samplingintv = ceil(nrows/newsize); %This is the sampling interval for block computation
                half_exp_size = ceil((blk_size(1)-samplingintv)/2)  %Calculate the expansion size of the block so that the
              
                lscal = @(block_struct) (4.0/var(reshape(block_struct.data,1,(2*half_exp_size+samplingintv)^2)));
                    %calculation is perform
                lsim = blockproc(org_im,[samplingintv samplingintv],lscal,'PadMethod','Replicate','BorderSize',...
                    [half_exp_size half_exp_size],'TrimBorder',0,'PadMethod','Replicate', 'PadPartialBlocks',1);
                
                %Now, calculate the gradient magnitude, Convert into the phase value
                org_im = imfilter(org_im,hs,'same'); %Smooth out the image before computing the gradient
                gx=imfilter(org_im,hx,'same');
                gy=imfilter(org_im,hy,'same');
                gmag = gx.^2+gy.^2;

                gmagcal = @(block_struct) mean(reshape(block_struct.data,1,(2*half_exp_size+samplingintv)^2));
                    %calculation is perform
                gmagavg = blockproc(gmag,[samplingintv samplingintv],gmagcal,'PadMethod','Replicate','BorderSize',...
                    [half_exp_size half_exp_size],'TrimBorder',0,'PadMethod','Replicate', 'PadPartialBlocks',1);

                clear org_im;
                clear gx;
                clear gy;
                save(strcat(curpath,ls_file_name),'curdist'); %Save the phase histogram distribution
                disp('Done.');
            end
            load(ls_file_name);
            figure(1);
            plot(curdist(1,:),'r');
            hold on;
            plot(curdist(2,:),'g');
            hold on;
            plot(curdist(3,:),'b');
            hold off;
            drawnow;
            legend('Lumen','Gland','Stroma');
       end
    end
    close(h);

end
