function [filters]=updateEnergy(filenames,curpath)
%This function computes the responses of images in the dataset w.r.t the
%odd and even symmetric filters
%Inputs:
%   filenames the name of all cores
%   path: the place to save phase distributions
%Outputs:
%   dist: a 4x1 cell array containing distributions corresponding to 4 different
%   classes. Each cell is a 3x500 matrix where each row is a histogram
%   distribution of the phase value...
    
    addpath(curpath);
    nrows = 3072;
    ncols=3072;
    addpath(strcat(cd(cd('..')),'\support\Texton_generator'));
    %addpath('C:\Users\Nikon\Dropbox\current working code\Texton_generator');
    
    %Parameters for filter response generator
    texton_num=20;
    texton_calc_en=0;%Enable/Disable texton calculation
    nscales = 5; ninvscales = 6; ndir=8;
    h = waitbar(0,'Filter response calc process...');
    retrain=0;
    nsamples = size(filenames,1);
    for sampleIdx=1:nsamples
            tic;
            waitbar(sampleIdx/nsamples,h,'Progress...')          
            %cur_file_name = filenames{classidx,1}{sampleIdx,1};
            cur_file_name = filenames{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);

            %cur_file_name = filenames{classidx,1}{sampleIdx,1};
            cur_file_name = strcat(cur_file_name(1:end-4),'_small.tif');%Read the small size image only
            fr_file_name = strcat(label_name,'_lm_fr.mat');
            disp(['Calculating LM filter response: ' label_name ' ...']);
            if ((~exist(fr_file_name,'file'))||(retrain==1))
                org_im = imread(cur_file_name); %Load the original data              
                
                im = imresize(org_im,[nrows ncols]);
                clear org_im;  

                %Generate the filter response and filter's psfs.
                [texton,texton_map,filters,f_res]=texton_compute(im,texton_num,'lm',0,'kmean',texton_calc_en);
                
                %nscales = 5; ninvscales = 6; ndir=8;
                %[fe,fo,emag,emap,edir]=computeOrientedEnergyLM(f_res,nscales,ninvscales,ndir);
                
                %Compute the maximum response feature
                disp('Calulating maximum response....');
                max_response_images = zeros(nrows,ncols,2*nscales+2*ninvscales);
                L = zeros(nrows,ncols);
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
                mr_data=zeros(nrows*ncols,2*nscales+2*ninvscales); %Each data is a row for segmentation
                for bandIdx=1:2*nscales+2*ninvscales
                    mr_data(:,bandIdx)=reshape(max_response_images(:,:,bandIdx),nrows*ncols,1);
                end
                L = L+sum(f_res(:,:,2*nscales*ndir+1:2*nscales*ndir+2*ninvscales).^2,3);
                L = sqrt(L);
                L = 1./L;
                save(strcat(curpath,fr_file_name),'mr_data','im','L','-v6'); %Save the phase histogram distribution - use v6 to save without compression. it will speed up the load speed.
                clear f_res;
                 %End of the code section for the paper...           
                 timeElapse = toc;
                 disp(['Calculation time: ' num2str(timeElapse) '(s)']);

            else
                %load(fr_file_name,'mr_data','L');
            end
    end
    close(h);

end