function computeTextonHistForWholeCore(filenames,labelpath,textonpath)
    %This function compute the 
    %This function computes the size of different glands in the biopsies
    %core...The glandsize as well as the distribution of glandsize will be
    %plotted...This code works on the 2048 x 2048 image. Otherwise, the
    %image has to be rescaled...
    %Compute the parameters of the biopsies using the groundtruth label...
    %To describe the stroma, we use the histogram of texton.
    addpath(textonpath);
    ntextons = 120;
    feat = zeros(0,ntextons);
    dataclass=zeros(0,1);
    filename = cell(0,1);
    color_arr='rgbyck';
    ngrades = size(filenames,1);
    for classidx=1:ngrades %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            tx_file_name = strcat(label_name,'_texton.mat'); %Load the texton map for current image
            lbl_file_name = strcat(label_name,'_resized.mat'); %Load the texton map for current image
                  
            disp(['Adding ' label_name ' ...']);
            if (~exist(strcat(textonpath,'diagtexton.mat'),'file'))                
                load(strcat(labelpath,lbl_file_name));         
                load(strcat(textonpath,tx_file_name));
                nonlumenidx=find(lblim~=0);
                textonhist=hist(new_text_map(nonlumenidx),[1:ntextons])/length(nonlumenidx);
                %Save the feature vector and the label...
                feat(end+1,:) = textonhist;
                dataclass(end+1,:)= classidx;
                filename{end+1,:}=label_name;
                figure(1);
                hold on;
                plot(feat(end,:),color_arr(mod(classidx,6)+1));
                %imagesc(feat);
            
            end
       end
    end
    
    save('diagtexton.mat','feat','dataclass','filename');
end