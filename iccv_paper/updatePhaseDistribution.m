function dist=updatePhaseDistribution(filenames,curpath,labelpath)
%This function computes the phase distributions for different classes
%Inputs:
%   filenames, glandnames: the name of all cores
%   path: the place to save phase distributions
%Outputs:
%   dist: a 4x1 cell array containing distributions corresponding to 4 different
%   classes. Each cell is a 3x500 matrix where each row is a histogram
%   distribution of the phase value...
    nbins = 200;
    dist = cell(4,1);
    for classidx=1:4
        dist{classidx,1}=zeros(3,nbins);
    end
    %part1: compute the distribution for each image and save it
    addpath(curpath);
    h = waitbar(0,'Calculating phase distribution progress...');
    histbin = linspace(0,65536,nbins);
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            curdist = zeros(3,nbins);
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            dist_file_name = strcat(label_name,'_phase_dist.mat');
            label_file_name = strcat(label_name,'.mat');
            if (~exist(dist_file_name,'file'))
                disp(['Calculating phase histogram: ' label_name ' ...']);
                org_im = imread(cur_file_name); %Load the original data
                load(strcat(labelpath,label_file_name));
                %Resize the phase and the label to new size
           %     org_im =cast(imresize(org_im,[2048 2048],'nearest'),'single');
           %     lblim = imresize(lblim,[2048 2048],'nearest');
                
                %Go through each class and compute the distributions
                for classIdx=1:3
                    pixidx = find(lblim==(classIdx-1));
                    histvect=hist(org_im(pixidx),histbin)/length(pixidx);
                    curdist(classIdx,:)=histvect(:)';%Assign to the histogram matrix
                end
                save(strcat(curpath,dist_file_name),'curdist'); %Save the phase histogram distribution
                disp('Done.');
            
            end
            load(dist_file_name);
            dist{classidx,1}=dist{classidx,1}+curdist; %Update the class distribution
            figure(1);
            plot(curdist(1,:),'r');
            hold on;
            plot(curdist(2,:),'g');
            hold on;
            plot(curdist(3,:),'b');
            hold on;
            drawnow;
            legend('Lumen','Gland','Stroma');
            
       end
       dist{classidx,1}=dist{classidx,1}/nsamples;
    end
    close(h);

end