function demo
    clc;
    clear all;
    close all;
    datapath = 'F:\QPI_data_fullset\';
    textpath = strcat(datapath,'textdir\');
    addpath(cd(cd('..')));
    svm_dir=strcat(datapath,'svm\ws20\');
    label_dir = strcat(datapath,'\metric&label\');
    load('texton_data.mat'); %Load the texton data
    ntextons=size(texton,1); %Get number of textons
    addpath(textpath); %Add the texton path      
    addpath(strcat(cd(cd('..\..')),'\Texton_generator'));
    addpath(strcat(cd(cd('..\..')),'\Bilateral filtering'));
    addpath(strcat(cd(cd('..')),'\trainedclassifier'));
    addpath(strcat(cd(cd('..')),'\support'));
    load 'svmstruct_ws40.mat';
   
    %Look for all file names
    filenames =  dir(strcat(datapath,'*.tif'));
    nfiles = length(filenames);
    filelist = cell(0,1);
    for filenameidx=1:nfiles
        filelist{end+1,1} =strcat(datapath,filenames(filenameidx,1).name);
    end
    
    h = waitbar(0,'Filter response calc process...');
    for fileIdx = 1:nfiles
        tic;
        cur_file_name = filelist{fileIdx,1};
        waitbar(fileIdx/nfiles,h,'Progress...')
        dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
        slash_pos = strfind(cur_file_name,'\');
        label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
        org_im = imread(cur_file_name); %Load the original data
        nrows = 2048;
        ncols = 2048;    
        npixels = nrows * ncols;
        radius = 40;

        org_im = imresize(org_im,[nrows ncols]);
        figure(1);
        imagesc(cast(org_im,'single')/65536*4.0-0.5);
        colorbar;
        drawnow;
        colormap jet;
        title('Original phase image');
        axis off;

        %First, calculate the filter response
        fr_file_name = strcat(label_name,'_lm_fr.mat');
        disp(['Calculating LM filter response: ' label_name ' ...']);
        if (~exist(fr_file_name,'file'))
            nscales = 5; ninvscales = 5; ndir=8;
            [texton,texton_map,filters,f_res]=texton_compute(org_im,20,'lm',0,'kmean',0);
            [fe,fo,emag,emap,edir]=computeOrientedEnergyLM(f_res,nscales,ninvscales,ndir);
            %Copute the maximum response
             disp('Calulating maximum response over the direction....');
             max_response_images = zeros(nrows,ncols,2*nscales+2*ninvscales);
             for scaleIdx=1:nscales
                   cur_even_data_block = f_res(:,:,(scaleIdx-1)*ndir+1:scaleIdx*ndir);
                   max_response_images(:,:,scaleIdx) = max(cur_even_data_block,[],3);
                   cur_odd_data_block = f_res(:,:,nscales*ndir+(scaleIdx-1)*ndir+1:nscales*ndir+...
                           scaleIdx*ndir);
                   max_response_images(:,:,nscales+scaleIdx) = max(cur_odd_data_block,[],3);
             end
             max_response_images(:,:,2*nscales+1:end)=f_res(:,:,2*nscales*ndir+1:2*nscales*ndir+2*ninvscales);
             mr_data=zeros(2048^2,2*nscales+2*ninvscales); %Each data is a row for segmentation
             for bandIdx=1:2*nscales+2*ninvscales
                  mr_data(:,bandIdx)=reshape(max_response_images(:,:,bandIdx),2048^2,1);
             end    
             save(strcat(textpath,fr_file_name),'fe','fo','emag','emap','edir','f_res','mr_data','-v7.3'); %Save the phase histogram distribution       
             clear max_response_images;
             clear f_res;
             clear fe;
             clear fo;
             clear emag;
             clear emap;
             clear edir;
        else
             load(fr_file_name,'mr_data');
        end
        figure(2);
        imagesc(reshape(mr_data(:,3),[nrows ncols]));
        axis off;
        drawnow;
        title('Maximum response at scale 3')

        %Next, update texton index for each pixel
        disp(['Performing vector quantization...']);
        tx_file_name = strcat(label_name,'_texton.mat'); %Name of the texton file
        ntextons=size(fin_texton,1);

        if (~exist(strcat(textpath,tx_file_name)))
            dist_map = 1e+20*ones(npixels,1); %This matrix store the distance from current filter response to all the new texton
            new_text_map=ones(nrows,ncols);
            for newtextonIdx=1:ntextons
                cur_dist = sum((mr_data-repmat(fin_texton(newtextonIdx,:),[npixels 1])).^2,2); %Compute the distance to the current label                        )).^2,2); %Computing distance from the set of current texton to the common set
                bestTextonIdx = find(cur_dist<dist_map);
                new_text_map(bestTextonIdx)=newtextonIdx;
                dist_map( bestTextonIdx)=cur_dist(bestTextonIdx);
                figure(3);
                imagesc(new_text_map);
                colorbar
                axis off;
                drawnow
                title('New indexing map');
            end  
            save(strcat(textpath,tx_file_name),'new_text_map','dist_map'); %Save the texton indexing image
            clear dist_map;
        else
            load(strcat(textpath,tx_file_name),'new_text_map','dist_map'); %Save the texton indexing image
        end
        figure(3);
        imagesc(new_text_map);
        colorbar
        axis off;
        drawnow;
        title('Texton indexing map');

        %Third step, computing the histogram of textons....
        disp(['Calculating histogram of textons....']);
        texton_hist_file_name =  strcat(label_name,'_texton_hist_',num2str(radius),'.mat');
       
        pixelidx=[1:npixels];
        pixelidx=pixelidx(:); %Line up current coordinates of the pixels
        x_offset = [-(radius/2) (radius/2) (radius/2) -(radius/2)];
        y_offset = [-(radius/2) -(radius/2) (radius/2) (radius/2)];
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
        if (~exist(strcat(textpath,texton_hist_file_name)))
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
            clear intIm;
            save(strcat(textpath,texton_hist_file_name),'histim','-v7.3');      
        else
            load(strcat(textpath,texton_hist_file_name),'histim');
        end
        clear pixelidx;
        clear nneighbors;
        clear neigh_coord;
        clear x_coord;
        clear y_coord;
        clear neigh_x;
        clear neigh_y;

        %Next, perform the segmentation between gland and stroma
        lblim = zeros(nrows,ncols);
        ds_mean_gland=mean(mean(org_im(1:30,1:30))); %Find the phase of the background region
        candidate_idx=find((org_im>1.08*ds_mean_gland)); %Make sure what we have is good for training
        lblim(candidate_idx)=1;
        se = strel('disk',6);
        lblim = imclose(lblim,se);


        %Compute the results on 512x512 images
        %Resize the image to the new size
        nrows = 512;
        ncols = 512; 
        npixels = nrows*ncols;
        lblim = imresize(lblim,[nrows ncols],'nearest');
        %Downsample the histogram of texton to the newsize
        
        svm_file_name = strcat(label_name,'_svm_200k_ws40.mat');
        if (~exist(strcat(svm_dir,svm_file_name)))
            fim = (-3)*ones(nrows,ncols);
            %Load the value of a classifier to compute
            shiftval = svmtruct.ScaleData.shift;
            scaleval = svmtruct.ScaleData.scaleFactor;
            supvect = svmtruct.SupportVectors;%Get the support vectors
            alpha = svmtruct.Alpha; %Note that alpha is positive for the 1st group and -1 for the second group
            bias = svmtruct.Bias;
            kerfunc = svmtruct.KernelFunction;
            kerfuncargs = svmtruct.KernelFunctionArgs;
            newhistim  = zeros(nrows,ncols,ntextons);
            for textonidx = 1:ntextons
                newhistim(:,:,textonidx) = imresize(histim(:,:,textonidx),[nrows ncols],'nearest');
            end
            histim = newhistim;

            %Line up to convert from 3D data to 2D data
            testdata = zeros(nrows*ncols,ntextons);
            for bandidx=1:ntextons
                curband = histim(:,:,bandidx);
                curvect = curband(:);
                testdata(:,bandidx)=curvect(:);
            end
          
            %Now, justclassify gland vs stroma, not go with lumen
            tempidx = find(lblim~=0);
            valset = testdata(tempidx,:);

            nsampleperbatch =4096;
            nbatch=ceil(size(valset,1)/nsampleperbatch);
            f=zeros(size(valset,1),1); %This is the value of the classification function
            for batchidx=1:nbatch
                 disp(['Batch index: ', num2str(batchidx), '/' num2str(nbatch)]);
                 curevalset =  valset((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,size(valset,1)),:);                       
                 curevalset = bsxfun(@plus,curevalset,shiftval);%Data normalization
                 curevalset = bsxfun(@times,curevalset,scaleval);
                 temp_out=-(kerfunc(supvect,curevalset,kerfuncargs{:})'*alpha(:) + bias(:));
                 f((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,size(valset,1)),:)=...
                    temp_out;
            end
            fim(tempidx)=f;
            newlblim = zeros(size(lblim));
            thresh = -0.03;
            candidateidx = find(lblim~=0);
            candidatelabel = (-1)*ones(size(candidateidx));
            glandidx = find(fim(candidateidx)>thresh);
            candidatelabel(glandidx)=1;
            newlblim(candidateidx)=candidatelabel;
            lumenidx = find(newlblim==0);
            newlblim(lumenidx)=-1;
            %Apply image dilation
            se = strel('disk',1);
            idxim = im2bw(newlblim,0);
            idxim = imclose(idxim,se);
            newlblim = cast(idxim,'single')*2.0-1.0;
            newlblim(lumenidx)=0;
            [newlblim] = procseg(newlblim); %Another processing step to remove sall are
            save(strcat(svm_dir,svm_file_name),'fim','f','newlblim');
        else
            load(strcat(svm_dir,svm_file_name),'fim','f','newlblim');
        end 
        clear histim;
        figure(4);
        imagesc(newlblim);
        colorbar
        axis off;
        clear fim;
        clear f;
        
        
        %Compute the results on 1024x1024 images
        %Resize the image to the new size
        nrows = 1024;
        ncols = 1024; 
        npixels = nrows*ncols;
        lblim = imresize(lblim,[nrows ncols],'nearest');
        %Downsample the histogram of texton to the newsize
        load 'svmstruct_ws40.mat'
        load(strcat(textpath,texton_hist_file_name),'histim');
        disp('Computing functional distance on 1024x1024 image')
        svm_file_name = strcat(label_name,'_svm_200k_ws40.mat');
        if (~exist(strcat(svm_dir,svm_file_name)))
            fim_1024 = (-3)*ones(nrows,ncols);
            %Load the value of a classifier to compute
            shiftval = svmtruct.ScaleData.shift;
            scaleval = svmtruct.ScaleData.scaleFactor;
            supvect = svmtruct.SupportVectors;%Get the support vectors
            alpha = svmtruct.Alpha; %Note that alpha is positive for the 1st group and -1 for the second group
            bias = svmtruct.Bias;
            kerfunc = svmtruct.KernelFunction;
            kerfuncargs = svmtruct.KernelFunctionArgs;
            newhistim  = zeros(nrows,ncols,ntextons);
            for textonidx = 1:ntextons
                newhistim(:,:,textonidx) = imresize(histim(:,:,textonidx),[nrows ncols],'nearest');
            end
            histim = newhistim;

            %Line up to convert from 3D data to 2D data
            testdata = zeros(nrows*ncols,ntextons);
            for bandidx=1:ntextons
                curband = histim(:,:,bandidx);
                curvect = curband(:);
                testdata(:,bandidx)=curvect(:);
            end
          
            %Now, justclassify gland vs stroma, not go with lumen
            tempidx = find(lblim~=0);
            valset = testdata(tempidx,:);

            nsampleperbatch =4096;
            nbatch=ceil(size(valset,1)/nsampleperbatch);
            f=zeros(size(valset,1),1); %This is the value of the classification function
            for batchidx=1:nbatch
                 disp(['Batch index: ', num2str(batchidx), '/' num2str(nbatch)]);
                 curevalset =  valset((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,size(valset,1)),:);                       
                 curevalset = bsxfun(@plus,curevalset,shiftval);%Data normalization
                 curevalset = bsxfun(@times,curevalset,scaleval);
                 temp_out=-(kerfunc(supvect,curevalset,kerfuncargs{:})'*alpha(:) + bias(:));
                 f((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,size(valset,1)),:)=...
                    temp_out;
            end
            fim_1024(tempidx)=f;
            newlblim = zeros(size(lblim));
            thresh = -0.03;
            candidateidx = find(lblim~=0);
            candidatelabel = (-1)*ones(size(candidateidx));
            glandidx = find(fim_1024(candidateidx)>thresh);
            candidatelabel(glandidx)=1;
            newlblim(candidateidx)=candidatelabel;
            lumenidx = find(newlblim==0);
            newlblim(lumenidx)=-1;
            %Apply image dilation
            se = strel('disk',1);
            idxim = im2bw(newlblim,0);
            idxim = imclose(idxim,se);
            newlblim = cast(idxim,'single')*2.0-1.0;
            newlblim(lumenidx)=0;
            [newlblim] = procseg(newlblim); %Another processing step to remove sall are
            save(strcat(svm_dir,svm_file_name),'fim_1024','f','newlblim','-append');
        else
            load(strcat(svm_dir,svm_file_name),'fim_1024','f','newlblim');
        end 
        clear histim;
        figure(4);
        imagesc(newlblim);
        colorbar
        axis off;
        clear fim;
        clear f;
        figure(3);
        imagesc(fim_1024);
        
        %Do bilateral filtering then save the resuls of bilateral filtering
        %on 512 x 512 image
        load(strcat(svm_dir,svm_file_name),'fim','f','newlblim','fimbifilter','newlblimbf');
        if (~exist('fimbifilter'))
              idx = find(fim~=(-3));
              fimbifilter = bifilter(fim,10,0.5*std(fim(idx)));
              %Save the bilateral filterin results
            save(strcat(svm_dir,svm_file_name),'fimbifilter','-append');
        end
        
        figure(5);
        imagesc([fim fimbifilter]);
        title('Original - Bilateral results...');
        newlblim = zeros(size(lblim));
      
        
        if (~exist('newlblimbf'))
            %Create a label map for stroma...Set it to a relatively
            %high threshold so that all stroma is detected at the price
            %of more error to gland...
            thresh =-0.029;
            candidateidx = find(lblim~=0);
            candidatelabel = (-1)*ones(size(candidateidx));
            glandidx = find(fimbifilter(candidateidx)>thresh);
            candidatelabel(glandidx)=1; %Assign gland pixels
            newlblim(candidateidx)=candidatelabel;
            lumenidx = find(newlblim==0); %Assign lumen pixels
            newlblim(lumenidx)=-1;
            %Clean the gland map a little bit
            newlblim = im2bw(newlblim,0);
            %Remove very small pores...
            newlbliminv = 1-newlblim;
            minporesize = 5;
            newlbliminv = bwareaopen(newlbliminv,minporesize,8); %Remove very small pores less...
            glandnoholes = 1-newlbliminv;
            %Slightly connect regions that are separated
            se = strel('disk',1);
            glandnoholes = imclose(glandnoholes,se);
            glandnoholes = bwareaopen(glandnoholes,5,8);
            newlblim = cast(glandnoholes,'single')*2.0-1.0;
            newlblim(lumenidx)=0;
            newlblimbf = newlblim;
            save(strcat(svm_dir,svm_file_name),'newlblimbf','-append');
        else
                 load(strcat(svm_dir,svm_file_name),'newlblim','newlblimbf');
        end
        figure(4);
        imagesc([newlblimbf]);
        axis on;
        title('label with BF filtering results');
        clear newlblimbf;
        clear fimbifilter; %Prepare for the next image
        %-------------------------------------------------------------------
        %Do bilateral filtering then save the resuls of bilateral filtering
        %on 1024 x 1024 image
        load(strcat(svm_dir,svm_file_name),'fim_1024','fimbifilter_1024','newlblimbf_1024');
        %if (~exist('fimbifilter_1024'))
              idx = find(fim_1024~=(-3));
              fimbifilter_1024 = bifilter(fim_1024,10,0.5*std(fim_1024(idx)));
              %Save the bilateral filterin results
            save(strcat(svm_dir,svm_file_name),'fimbifilter_1024','-append');
        %end
        
        figure(5);
        imagesc([fim_1024 fimbifilter_1024]);
        title('Original - Bilateral results...');
        newlblim_1024 = zeros(size(lblim));
        
        
        %if (~exist('newlblimbf_1024'))
            %Create a label map for stroma...Set it to a relatively
            %high threshold so that all stroma is detected at the price
            %of more error to gland...
            thresh =-0.029;
            candidateidx = find(lblim~=0);
            candidatelabel = (-1)*ones(size(candidateidx));
            glandidx = find(fimbifilter_1024(candidateidx)>thresh);
            candidatelabel(glandidx)=1; %Assign gland pixels
            newlblim_1024(candidateidx)=candidatelabel;
            lumenidx = find(newlblim_1024==0); %Assign lumen pixels
            newlblim_1024(lumenidx)=-1;
            %Clean the gland map a little bit
            newlblim_1024 = im2bw(newlblim_1024,0);
            %Remove very small pores...
            newlbliminv_1024 = 1-newlblim_1024;
            minporesize = 5;
            newlbliminv_1024 = bwareaopen(newlbliminv_1024,minporesize,8); %Remove very small pores less...
            glandnoholes = 1-newlbliminv_1024;
            %Slightly connect regions that are separated
            se = strel('disk',1);
            glandnoholes = imclose(glandnoholes,se);
            glandnoholes = bwareaopen(glandnoholes,5,8);
            newlblim_1024 = cast(glandnoholes,'single')*2.0-1.0;
            newlblim_1024(lumenidx)=0;
            newlblimbf_1024 = newlblim_1024;
            save(strcat(svm_dir,svm_file_name),'newlblimbf_1024','-append');
            %Save the label image as well as the bilaterial filtering in an
            %image file
            writeTIFF(fimbifilter_1024,strcat(label_dir,label_name,'_met.tif'));
            writeTIFF(newlblimbf_1024,strcat(label_dir,label_name,'_lbl.tif'));
         
        %else
             load(strcat(svm_dir,svm_file_name),'newlblim_1024','newlblimbf');
        %end
        figure(4);
        imagesc([newlblimbf_1024]);
        axis on;
        title('label with BF filtering results')
        clear newlblimbf_1024;
        clear fimbifilter_1024;
        avgTime = toc;
        disp(['Time spent for 1 sample: ' num2str(avgTime)]);
    end
    close(h);
    
end






