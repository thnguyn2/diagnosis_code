function [texton]=updateTexton(filenames,tifflist,labelpath,curpath)
%This function computes the texton of single images and a final texton for
%all the images. The map of corresponding texton is also included
%odd and even symmetric filters
%Inputs:
%   filenames the name of all cores
%   curpath: path to the directory that saves texton information
%   labelpath: path to the directory that saves the label information.

%Outputs:
%   dist: a 4x1 cell array containing distributions corresponding to 4 different
%   classes. Each cell is a 3x500 matrix where each row is a histogram
%   distribution of the phase value...
    
    addpath(curpath);
    nrows = 3072;
    ncols = 3072;
    %Parameters for filter response generator
    nscales = 5; ninvscales = 5; 
    nspectraperclass = 10000;
    nbands = 20;

    filter_response = zeros(0,nbands);
    label_vect = zeros(0,nbands);
 
%    texton_num_array=round(logspace(log10(5),log10(2000),30)); %This is the number of textons for training, we just train with different k and lets see the compactness
    texton_num_array = 50; %The best compromising point between speed and accuracy
    
%     texton_map_over_images = zeros(0,2*(nscales+ninvscales));
%     texton_label_over_images = zeros(0,2*(nscales+ninvscales));
%     prior_over_images = zeros(0,3);
%     hist_over_images = zeros(0,3);
    
    ngrades = size(filenames,1);
    retrain=0;
    %Now, compute the texton for each image
    if ((~exist('texton_data_before_kmeans.mat','file'))||(retrain==1))          
     for classidx=1:ngrades %Go through different classes
        nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
        for sampleIdx=1:nsamples
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            fr_file_name = strcat(label_name,'_lm_fr.mat');%Name of the filter response
            tx_file_name = strcat(label_name,'_texton.mat');
            label_file_name = strcat(label_name,'_resized.mat');
            disp(['Adding random samples of: ' label_name ' ...']);
            load(strcat(curpath,fr_file_name),'mr_data');
            %Load the label image
            load(strcat(labelpath,label_file_name));
            glandidx = find(lblim==1);
            stromaidx = find(lblim==2);
            lumenidx = find(lblim==0);
                
            nglandpix = length(glandidx); %Get number of gland and stroma pixs
            nstromapix = length(stromaidx);
            nlumenpix = length(lumenidx);
            selectedglandidx = randperm(nglandpix);
            selectedstromaidx = randperm(nstromapix);
            selectedlumenidx = randperm(nlumenpix);
            gland_subset=mr_data(...
            glandidx(selectedglandidx(1:nspectraperclass)),1:nbands);
            stroma_subset=mr_data(...
            stromaidx(selectedstromaidx(1:nspectraperclass)),1:nbands);
            lumen_subset=mr_data(...
            lumenidx(selectedlumenidx(1:nspectraperclass)),1:nbands);                
                          
            filter_response(end+1:end+nspectraperclass,:)=lumen_subset;
            filter_response(end+1:end+nspectraperclass,:)=gland_subset;
            filter_response(end+1:end+nspectraperclass,:)=stroma_subset;
            label_vect(end+1:end+nspectraperclass)=0;
            label_vect(end+1:end+nspectraperclass)=1;
            label_vect(end+1:end+nspectraperclass)=2;
                %save(strcat(curpath,tx_file_name),'gland_subset','stroma_subset','lumen_subset');
                
                
       end
     end
     save(strcat(curpath,'texton_data_before_kmeans.mat'),'filter_response','label_vect');
   end
     
      
   
       
    disp('Finding cluster centers...');
    retrain=0;
    computing_platform=1; %OpenCV-Mex functions
    nk = length(texton_num_array); %Number of k values
    compactness_arr = zeros(nk,1);
    %form filter response over images
    if ((~exist('texton_data_before_kmeans.mat','file'))||(retrain==1))
        load(strcat(curpath,'texton_data_before_kmeans.mat'));
        for kidx=1:nk
            tic
            texton_num = texton_num_array(kidx);
            disp(['Working on k= ' num2str(texton_num) ', ' num2str(round(nk+1-kidx)) ' values to go']);
            k_means_file_name = strcat(curpath,'kmeans_res_',num2str(texton_num),'_clusters.mat');
            if (computing_platform==0) %Matlab k-means, which is 20x-60x slower than OpenCVMex k-means
                opts = statset('Display','iter', 'MaxIter',3000);
                
                [newidxmap,texton,compactness]=kmeans(filter_response, texton_num,'Options',opts,'Replicates',1,'Start','sample');
                
            else
                [newidxmap, texton,compactness] = kmeans(filter_response,texton_num,'Attempts',10); %Compactness is the empirical error
               
            end
            compactness_arr(kidx)=compactness;
            figure(1);
            plot(texton_num_array(1:kidx),compactness_arr(1:kidx));
            title('Empirical vs k');
            save(k_means_file_name,'texton','newidxmap','compactness');
            timeElapse= toc
            disp(['Training time: ' num2str(timeElapse), ' (s)']);
        end
    end
  
    %Load the dictionary of texton and find the best matching index for
    %each pixel
    load(strcat(curpath,'kmeans_res_',num2str(texton_num_array),'_clusters.mat'));
     ntextons=size(texton,1);       
     retrain = 0;
     %Vector quantization on the dictionary of the overall textons. This
     %has to be done for all the tiff files that are available
     nsamples = size(tifflist,1);
     for sampleidx=1:nsamples
   %  for classidx=1:ngrades %Go through different classes
   %        nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
   %        for sampleIdx=1:nsamples
                   cur_file_name = tifflist{sampleidx,1};
                   %cur_file_name = filenames{classidx,1}{sampleIdx,1};
                   dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                   slash_pos = strfind(cur_file_name,'\');
                   label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                   tx_file_name = strcat(label_name,'_texton.mat'); %Load the texton map for current image
                   fr_file_name = strcat(label_name,'_lm_fr.mat');%Name of the filter response
                                    
                   disp(['Re-assigning index based on general textons: ' label_name ' ...']);
                   if ((~exist(strcat(curpath,tx_file_name)))||(retrain==1))
                       tic
                       load(strcat(curpath,fr_file_name),'mr_data'); %Load the current texton map and the maximum response data
                       npixels = size(mr_data,1);
                  
                       if (computing_platform==0)
                           dist_map = 1e+20*ones(npixels,1); %This matrix store the distance from current filter response to all the new texton
                           new_text_map=ones(2048,2048);
                           for newtextonIdx=1:ntextons
                               cur_dist = sum((mr_data(:,1:nbands)-repmat(texton(newtextonIdx,:),[npixels 1])).^2,2); %Compute the distance to the current label                        )).^2,2); %Computing distance from the set of current texton to the common set
                               bestTextonIdx = find(cur_dist<dist_map);
                               new_text_map(bestTextonIdx)=newtextonIdx;
                               dist_map( bestTextonIdx)=cur_dist(bestTextonIdx);                           
                           end
                       else
                           %Compute the label with openCV using 1-NN - 12x
                           %faster. 10 s/image
                           kn=KNearest(texton,[1:ntextons]);
                           new_text_map=kn.predict(mr_data(:,1:nbands));
                           new_text_map=reshape(new_text_map,[nrows ncols]);
                       end                   
                       timeElap = toc;
                       disp(['Vector quantization time: ' num2str(timeElap)]);
                   else
                       load(strcat(curpath,tx_file_name));
                   end
                   figure(1);
                   imagesc(new_text_map);
                   colorbar
                   drawnow
                   title('New indexing map');
                   save(strcat(curpath,tx_file_name),'new_text_map'); %Save the texton indexing image 
     %      end
     end

  
end


    