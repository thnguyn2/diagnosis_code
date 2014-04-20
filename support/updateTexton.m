function [texton]=updateTexton(filenames,labelpath,curpath)
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
    addpath('C:\Users\thnguyn2\Dropbox\current working code\Texton_generator');
    
    %Parameters for filter response generator
    texton_num_per_im=100;
    nscales = 5; ninvscales = 5; ndir=8;
    subsetsize = 50000; %Randomly select the pixel

    h = waitbar(0,'Calculating textons...');
    k = waitbar(0,'Step progress...');
    texton_map_over_images = zeros(0,2*(nscales+ninvscales));
    texton_label_over_images = zeros(0,2*(nscales+ninvscales));
    prior_over_images = zeros(0,3);
    hist_over_images = zeros(0,3);
    %% Computing textons for each images..... 
    nfilestotal = 0;
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            nfilestotal =nfilestotal+1; 
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            max_response_images = zeros(2048,2048,2*nscales+2*ninvscales);
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            fr_file_name = strcat(label_name,'_lm_fr.mat');%Name of the filter response
            tx_file_name = strcat(label_name,'_texton.mat');
            if (~exist(fr_file_name,'file'))
                disp(['Calculating textons for: ' label_name ' ...']);

                load(strcat(curpath,fr_file_name),'f_res');
                %Compute the maximum response over non-symmetric filers
                L=zeros(2048,2048);
                disp('Calulating maximum response....');
                for scaleIdx=1:nscales
                    cur_even_data_block = f_res(:,:,(scaleIdx-1)*ndir+1:scaleIdx*ndir);
                    max_response_images(:,:,scaleIdx) = max(cur_even_data_block,[],3);
                    L=L+sum(cur_even_data_block.^2,3);
                    cur_odd_data_block = f_res(:,:,nscales*ndir+(scaleIdx-1)*ndir+1:nscales*ndir+...
                        scaleIdx*ndir);
                    L=L+sum(cur_odd_data_block.^2,3);
                    max_response_images(:,:,nscales+scaleIdx) = max(cur_odd_data_block,[],3);
                end
                max_response_images(:,:,2*nscales+1:end)=f_res(:,:,2*nscales*ndir+1:2*nscales*ndir+2*ninvscales);
                
                %There is no need to perform normalization here since our
                %measurement are quantitative
%                L = L+sum(f_res(:,:,2*nscales*ndir+1:2*nscales*ndir+2*ninvscales).^2,3);
%                L = sqrt(L);
%                L = 1./L;

                mr_data=zeros(2048^2,2*nscales+2*ninvscales); %Each data is a row for segmentation
                for bandIdx=1:2*nscales+2*ninvscales
                    %No need to do normalization here...
                    %max_response_images(:,:,bandIdx)=max_response_images(:,:,bandIdx).*L;
                    mr_data(:,bandIdx)=reshape(max_response_images(:,:,bandIdx),2048^2,1);
                end
                save(strcat(curpath,fr_file_name),'mr_data','-append');
            end
%       load(strcat(curpath,tx_file_name),'text_map','texton_map','mr_data','idxmap');
       end
    end
   
    %Now, compute the texton for each image
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
            fr_file_name = strcat(label_name,'_lm_fr.mat');%Name of the filter response
            tx_file_name = strcat(label_name,'_texton.mat');
            label_file_name = strcat(label_name,'.mat');
            tx_file_name = strcat(label_name,'_texton.mat');
            if (~exist(tx_file_name,'file'))
                disp(['Calculating textons for: ' label_name ' ...']);
                load(strcat(curpath,fr_file_name),'mr_data');
                %Load the label image
                load(strcat(labelpath,label_file_name));
                lblim = imresize(lblim,[2048 2048],'nearest');
                glandidx = find(lblim==1);
                stromaidx = find(lblim==2);
                lumenidx = find(lblim==0);
                disp('Finding cluster centers....');
                nglandpix = length(glandidx); %Get number of gland and stroma pixs
                nstromapix = length(stromaidx);
                nlumenpix = length(lumenidx);
                selectedglandidx = randperm(nglandpix);
                selectedstromaidx = randperm(nstromapix);
                selectedlumenidx = randperm(nlumenpix);
                gland_subset=mr_data(...
                    glandidx(selectedglandidx(1:subsetsize)),:);
                stroma_subset=mr_data(...
                    stromaidx(selectedstromaidx(1:subsetsize)),:);
                lumen_subset=mr_data(...
                    lumenidx(selectedlumenidx(1:subsetsize)),:);
                save(strcat(curpath,tx_file_name),'gland_subset','stroma_subset','lumen_subset');
                %Save subsamples of the data

                
            end
        end
     end
      
    
       
    close(k);
      
    k = waitbar(0,'Step progress...');

    %Next, perform a secondary grouping over the texton
    disp('Grouping of texton...');
    
    nspectraperclass = 750;
    nbands = 20;
    filter_response = zeros(0,nbands);
    label_vect = zeros(0,nbands);
    %form filter response over images
    if (~exist('texton_data.mat','file'))
         for classidx=1:4 %Go through different classes
             nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
             for sampleIdx=1:nsamples
                 waitbar(sampleIdx/nsamples,h,'Progress...');
                 cur_file_name = filenames{classidx,1}{sampleIdx,1};
                 %Check to see if the label file exist
                 dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                 slash_pos = strfind(cur_file_name,'\');
                 label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                 tx_file_name = strcat(label_name,'_texton.mat');
 
                 disp(['Loading textons information for: ' label_name ' ...']);
                 load(strcat(curpath,tx_file_name),'gland_subset','stroma_subset','lumen_subset');
                 filter_response(end+1:end+nspectraperclass,:)=gland_subset(1:nspectraperclass,:);
                 filter_response(end+1:end+nspectraperclass,:)=stroma_subset(1:nspectraperclass,:);
                 label_vect(end+1:nspectraperclass)=1;
                 label_vect(end+1:nspectraperclass)=2;
            end
         end
        texton_num = 300;
        opts = statset('Display','iter', 'MaxIter',3000);
        [newidxmap,texton]=kmeans(filter_response, texton_num,'Options',opts,'Replicates',1,'Start','sample');
        save('texton_data.mat','texton','newidxmap','filter_response','label_vect');
    end
    [texton]=optimize_kmeans();
    ntextons=size(texton,1);
    close all;
    
    %Vector quantization on the dictionary of the overall textons
    for classidx=1:4 %Go through different classes
          nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
          for sampleIdx=1:nsamples
                  waitbar(sampleIdx/nsamples,h,'Calculating new indexing images...')
                  cur_file_name = filenames{classidx,1}{sampleIdx,1};
                  dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                  slash_pos = strfind(cur_file_name,'\');
                  label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                  tx_file_name = strcat(label_name,'_texton.mat'); %Load the texton map for current image
                  fr_file_name = strcat(label_name,'_lm_fr.mat');%Name of the filter response
                  load(strcat(curpath,fr_file_name),'mr_data'); %Load the current texton map and the maximum response data
                  npixels = size(mr_data,1);
                 
                  disp(['Re-assigning index based on general textons: ' label_name ' ...']);
                  dist_map = 1e+20*ones(npixels,1); %This matrix store the distance from current filter response to all the new texton
                  new_text_map=ones(2048,2048);
                  for newtextonIdx=1:ntextons
                      cur_dist = sum((mr_data-repmat(texton(newtextonIdx,:),[npixels 1])).^2,2); %Compute the distance to the current label                        )).^2,2); %Computing distance from the set of current texton to the common set
                      bestTextonIdx = find(cur_dist<dist_map);
                      new_text_map(bestTextonIdx)=newtextonIdx;
                      dist_map( bestTextonIdx)=cur_dist(bestTextonIdx);
                      figure(1);
                      imagesc(new_text_map);
                      colorbar
                      drawnow
                      title('New indexing map');
                  end           
                  %Now, compute the distribution corresponding to each
                  %texton. 
                 lbl_file_name = strcat(label_name,'.mat');
                 load(strcat(labelpath,lbl_file_name));
                 %Downsample the images to match the size
                 lblimnew = imresize(lblim,[2048 2048],'nearest');
                 %Find the number of texton that in each class and perform
                 %pruning...
                 
                 %Compute the prior distribution
                 nlumen = length(find(lblimnew==0));
                 ngland = length(find(lblimnew==1));
                 nstroma= length(find(lblimnew==2));
                 ovr_prior = [nlumen/npixels ngland/npixels nstroma/npixels];
                 ovr_posterior = zeros(ntextons,3);
                 for textonIdx = 1:ntextons
                      pixidx = find(new_text_map==textonIdx);
                      histo=hist(lblimnew(pixidx),[0 1 2])./[nlumen ngland nstroma];
                      ovr_posterior(textonIdx,:)=histo;
                 end
                 figure(2);
                 imagesc(ovr_posterior);
                 colorbar
                 drawnow;
                 title('Overall posterior')
                 figure(3);
                 imagesc(lblimnew);
                 save(strcat(curpath,tx_file_name),'new_text_map','ovr_prior','ovr_posterior','-append'); %Save the texton indexing image

          end
    end
   
end

function [fin_texton,fin_datacount]=optimize_kmeans()
    load('texton_data.mat','texton','newidxmap','filter_response','label_vect');
    ntextons = size(texton,1);
    nbands = size(texton,2);
    %Count how many data associate with each textons
    datacountpertexton = hist(newidxmap,[1:ntextons]);
    ndatasamples = size(filter_response,1);
    
    %Round 1, elliminates the textons that has the number of counts less then
    %the average
    textonnumthresh = round(ndatasamples/ntextons);
    eli_texton_idx = find(datacountpertexton<textonnumthresh);%Get the index of eliminated threshold
    ret_texton_idx = find(datacountpertexton>=textonnumthresh);
    %Compute the anfinity matrix between two textons
    d=affinity_mat(texton);
    %Maximize the elements in the diagonal
    d = d+eye(ntextons)*max(d(:));
    
    ret_texton = texton(ret_texton_idx,:);%This is the set of retained textons
    ret_datacount =  datacountpertexton(ret_texton_idx);
    
    %Go through all textons and merge the eliminated one to another one
    %that is closest to it
    neli = length(eli_texton_idx);
    
    not_consider_list = zeros(0,1); %This list contains the testons that has already been grouped and will not be considered
    further_consider_texton = zeros(0,nbands);
    for eliidx=1:neli %Go through
        %Find the closest textons to merge
        textoneliidx=eli_texton_idx(eliidx);
        cur_dist = d(textoneliidx,:); %Get the distance from the current textons to all of its neighbors
        %Mask out all the textons in the eliminated list and will not
        %consider it for grouping by set the distance to max value
        cur_dist(eli_texton_idx)=1e+6; %Set the distance corresponds to of eliminated guys to a dummy large number
        [closest_dist,bestmergeidx]=min(cur_dist);
        %Update the texton's coordinates
        eli_texton=texton(textoneliidx,:);
        merge_texton =texton(bestmergeidx,:); 
        cur_texton = (merge_texton* datacountpertexton(bestmergeidx)+eli_texton* datacountpertexton(textoneliidx))/...
            (datacountpertexton(bestmergeidx) + datacountpertexton(textoneliidx));
        merge_texton_idx_in_ret_set = find(ret_texton_idx==bestmergeidx);
        ret_texton(merge_texton_idx_in_ret_set,:)=cur_texton;
        ret_datacount(merge_texton_idx_in_ret_set)=ret_datacount(merge_texton_idx_in_ret_set)+datacountpertexton(textoneliidx);
        textonidxofpixelinnewdataset=find(newidxmap==textoneliidx);%Reupdate the index to the new in the retained dataset
        newidxmap(textonidxofpixelinnewdataset)=bestmergeidx;
        figure(1);
        plot(eli_texton,'r');
        hold on;
        plot(merge_texton,'g');
        hold on;
        plot(cur_texton,'b');
        hold off;
        figure(2);
        stairs(ret_datacount);
        title('Counts per texton');
        
    end
    ret_newidxmap = zeros(size(newidxmap));
    %Recompute the indexing map
    for rettexidx=1:length(ret_texton_idx)
        tempidx=find(newidxmap==ret_texton_idx(rettexidx));
        ret_newidxmap(tempidx)=rettexidx;
    end

    texton = ret_texton;
    datacount = ret_datacount;
    newidxmap = ret_newidxmap;
    ntextons = size(texton,1);
 
    %Compute the ratio of glands vs stroma for each textons...  
    
    %Round 2, group all textons that are very close together and have the
    %don't have very high specificity
    %Look for the closest cluster centers
    d=affinity_mat(texton);
    %Maximize the elements in the diagonal
    d = d+eye(ntextons)*max(d(:));
    allowdist = 0.4;
    mergecoupleidx =find(d<allowdist);
    
    [text1idx,text2idx]=ind2sub([ntextons,ntextons],mergecoupleidx);    
    %Just consider 1 couple once in which each elemenents is merged once to the
    %closest vectors
    
    %Now, form a list of groups containing similar textons
    text2idxunique = unique(text2idx);
    cluster_list =  cell(0,1);
    remaininglist = text2idxunique;
    deployed_list=zeros(0,1); %This list will contains all the textons that has been used in grouping to avoid the case that some
    %texton will appear in two groups
    while (1)
        cluster_list{end+1,1}=zeros(0,1); %Initialize a new group if still continue
        indexofclusterseed = remaininglist(1);
        cluster_list{end,1}(end+1,1)=indexofclusterseed;
        deployed_list(end+1,1)=indexofclusterseed;
        tempidx = find(text2idx==indexofclusterseed);%Look for the close guy in text1idx
        neighborcandidatelist = setdiff(text1idx(tempidx),deployed_list);
        nnbors = length(neighborcandidatelist);
        if (nnbors>0)
            cluster_list{end,1}(end+1:end+nnbors,1)=neighborcandidatelist;
            deployed_list(end+1:end+nnbors,1)=neighborcandidatelist;
        end
        %Now, continue with the items in the remaining set
        cur_cluster = cluster_list{end,1};
        remaininglist = setdiff(remaininglist,cur_cluster);
        if (isempty(remaininglist))
            break;
        end
    end
    
    nonupdatetexton =setdiff([1:ntextons],deployed_list);
    fin_texton = zeros(0,nbands);
    fin_datacount=zeros(0,1);
    fin_texton(end+1:end+length(nonupdatetexton),:)= texton(nonupdatetexton,:);
    fin_datacount(end+1:end+length(nonupdatetexton),1) = datacount(nonupdatetexton);
    ngroups=size(cluster_list,1);
    for groupidx=1:ngroups
        cur_text_list = cluster_list{groupidx,1};%Get a list of current textons
        cur_datacount = datacount(cur_text_list);
        cur_texton = texton(cur_text_list,:);
        new_texton = sum(cur_texton.*repmat(cur_datacount(:),[1 nbands]),1)/sum(cur_datacount);
        fin_texton(end+1,:)=new_texton;
        fin_datacount(end+1,1)=sum(cur_datacount);
    end
    
    ntextons=size(fin_texton,1);
    figure(1);
    hold off;
    %Display all textons and their histogram
    for textonidx=1:ntextons
        plot(fin_texton(textonidx,:));
        hold on;
    end
    figure(2)
    stairs(fin_datacount);
    title('Texton histogram')
    save('texton_data.mat','fin_texton','fin_datacount','-append');
    
end

function d = affinity_mat(texton)
    %Weight the contribution of each band by dividing the difference by
    %standard deviation at each band
    std_array = std(texton,0,1);
    ntextons = size(texton,1);
    nbands = size(texton,2);
    d = zeros(ntextons,ntextons);
    weighing_matrix = repmat(std_array,[ntextons 1]);
    for textonidx=1:ntextons
        distance = (repmat(texton(textonidx,:),[ntextons 1])-texton);
        d(textonidx,:)=(sum(abs(distance)./weighing_matrix,2)/nbands);
    end
end
 


    