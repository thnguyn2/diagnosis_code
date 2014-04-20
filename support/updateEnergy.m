function updateEnergy(filenames,curpath)
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
    addpath(strcat(cd(cd('..')),'\support\Texton_generator'));
    %addpath('C:\Users\Nikon\Dropbox\current working code\Texton_generator');
    
    %Parameters for filter response generator
    texton_num=20;
    texton_calc_en=0;%Enable/Disable texton calculation
    nscales = 5; ninvscales = 6; ndir=8;

    h = waitbar(0,'Filter response calc process...');
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            cur_file_name = strcat('H:\QPI data\',label_name,'.tif');
            fr_file_name = strcat(label_name,'_lm_fr.mat');
            %fr_file_name_fs = strcat(label_name,'_lm_fr_fs.mat'); %Full size image after filtering on 10000x10000 image
            disp(['Calculating LM filter response: ' label_name ' ...']);
            org_im = imread(cur_file_name); %Load the original data              
                
            im = imresize(org_im,[2048 2048]);
            clear org_im;  
            if (~exist(fr_file_name,'file'))
          
                %Generate the filter response and filter's psfs.
                [texton,texton_map,filters,f_res]=texton_compute(im,texton_num,'lm',0,'kmean',texton_calc_en);
                [fe,fo,emag,emap,edir]=computeOrientedEnergyLM(f_res,nscales,ninvscales,ndir);
                %save(strcat(curpath,fr_file_name),'fe','fo','emag','emap','edir','f_res','-v7.3'); %Save the phase histogram distribution
                clear f_res;
                clear fe;
                clear fo;
                clear emag;
                clear emap;
                disp('Done.');
            else
                load(fr_file_name,'fe','fo','emag','emap','edir');
            end
            
            %This section of the code is for the paper
             xl = 1000;xr = 1500;
             yt = 700;yb = 1200;
             figure(1);
             imagesc(cast(im(yt:yb,xl:xr),'single')*4.0/65536-0.5);
             colorbar;
             truesize;
             title('Phase image');
             axis off;
             figure(2);
             imagesc(fe(yt:yb,xl:xr,5)+fo(yt:yb,xl:xr,5));
             colorbar;
             truesize;
             title('Oriented energy');
             axis off;
             figure(3);
             imagesc((edir(yt:yb,xl:xr,5)-1)/8*180);
             colorbar;
             truesize;
             title('Max dir image');
             axis off;
             %End of the code section for the paper...
             
%             if (~exist(fr_file_name_fs,'file'))
%                 disp(['Calculating LM filter response [full size]: ' label_name ' ...']);
%                 org_im = imread(cur_file_name); %Load the original data              
%                 
%                 im = org_im;
%                 %Generate the filter response and filter's psfs.
%                 fs_enable = 1; %Enable the calculation of full-size images
%                 [texton,texton_map,filters,f_res]=texton_compute(im,texton_num,'lm',0,'kmean',texton_calc_en,fs_enable);
%                 [fe,fo,emag,emap,edir]=computeOrientedEnergyLM(f_res,nscales,ninvscales,ndir);
%                 save(strcat(curpath,fr_file_name_fs),'fe','fo','emag','emap','edir','f_res','-v7.3'); %Save the phase histogram distribution
%                 clear f_res;
%                 clear fe;
%                 clear fo;
%                 clear emag;
%                 clear emap;
%                 disp('Done.');
%             else
%                 load(fr_file_name_fs,'fe','fo','emag','emap','edir');
%             end
       end
    end
    close(h);

end