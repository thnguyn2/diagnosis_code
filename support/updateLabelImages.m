%Last update: April 19th
function updateLabelImages(filenames,glandnames,curpath)
    %This function computes the label images of all files in the dataset
    %and save it to the curpath.
    addpath(curpath);
    h = waitbar(0,'Calculating progress...');
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            cur_gland_file_name = glandnames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            label_file_name = strcat(label_name,'_resized.mat');
            disp(['Labeling: ' label_name ' ...']);
            if (~exist(label_file_name,'file'))
               
                org_im = imread(cur_file_name);
                gland_im = imread(cur_gland_file_name);
                lblim=findLabelImg(org_im,gland_im);
                %save(strcat(curpath,label_file_name),'lblim');               
            else
                load(strcat(curpath,label_file_name));
            end
            figure(1);
            imagesc(lblim);
        end
    end
    close(h);

end