function updateWeight(filenames,curpath,weighpath,cursize)
%This function computes the affinitiy matrix based on edge amplitude, edge direction and texton distribution difference
%at the required size (cursize)
%Input data:
%   filenames: a list of file names to compute the affinitiy matrix
%   curpath: path to the locations of features
%   cursize: size of the features (256/512/1024)
%   weighpath: path to the folder where the weights get stored.
%Output data:
%   [none]
%-------------------------------------
    addpath(curpath);
    addpath('C:\Users\thnguyn2\Dropbox\current working code\Texton_generator');
    nfilestotal = 0; %Number of files processed
    h = waitbar(0,'Calculating weighs for N-CUT...');
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
            
            %tx_file_name = strcat(label_name,'_texton.mat');
            
            %First, compute the weights based on the texture direction
            load(fr_file_name); %Load information about the filter response
           
       end
    end
    close(h);
end