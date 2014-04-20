function updateFilterResForNucleiCalc(filenames,curpath)
%This function computes the responses of images in the dataset w.r.t the
%odd and even symmetric filters. This section has finer filter responses at
%more angles
%Inputs:
%   filenames the name of all cores
%   path: the place to save phase distributions
%Outputs:
%   dist: a 4x1 cell array containing distributions corresponding to 4 different
%   classes. Each cell is a 3x500 matrix where each row is a histogram
%   distribution of the phase value...
    
    addpath(curpath);
    addpath('C:\Users\thnguyn2\Dropbox\current working code\Texton_generator');
    
    %Parameters for filter response generator
    texton_calc_en=0;%Enable/Disable texton calculation
    nscales = 1; ninvscales = 5; ndir=20;

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
            fr_file_name = strcat(label_name,'_lm_fr_fine_scale.mat');
            if (~exist(fr_file_name,'file'))
                disp(['Calculating LM filter response at finner scale: ' label_name ' ...']);
                org_im = imread(cur_file_name); %Load the original data
                
                
                im = imresize(org_im,[2048 2048]);
                clear org_im;  
                %Generate the filter response and filter's psfs.
                [texton,texton_map,filters,f_res]=texton_compute_for_nuclei(im,0,'lm',0,'kmean',texton_calc_en);
                [fe,fo,emag,emap,edir]=computeOrientedEnergyLM(f_res,nscales,ninvscales,ndir);
                save(strcat(curpath,fr_file_name),'fe','fo','emag','emap','edir','f_res','-v7.3'); %Save the phase histogram distribution
                disp('Done.');
            else
                load(fr_file_name,'fe','fo','emag','emap','edir');
            end
             figure(1)
             imagesc(emap(:,:,1).*edir(:,:,1));
             colorbar;
             drawnow;
       end
    end
   
  close(h);

end