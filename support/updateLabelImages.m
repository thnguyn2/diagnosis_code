%Last update: April 19th
function updateLabelImages(filenames,glandnames,curpath,corediagpath)
    %This function computes the label images of all files in the dataset
    %and save it to the curpath.
    addpath(curpath);
    h = waitbar(0,'Calculating progress...');
    ngrades = size(filenames,1);
    retrain = 0;
   
    for classidx=1:ngrades %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            cur_gland_file_name = glandnames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            label_file_name = strcat(label_name,'_label.tif');
            disp(['Labeling: ' label_name ' ...']);
            if ((~exist(label_file_name,'file'))|(retrain==1))
               
                org_im = imread(cur_file_name);
                gland_im = imread(cur_gland_file_name);
                lblim=findLabelImg(org_im,gland_im);
                core_lblim_name = strcat(corediagpath,label_name,'_roi_diag.tif');
                if (exist(core_lblim_name))
                    core_lblim = imread(core_lblim_name);
                    scrtidx = find((core_lblim==7.0)|(core_lblim==6.0));%index of scretion
                    if (~isempty(scrtidx))
                        lblim(scrtidx)=3.0;
                    end
                    bldidx = find(core_lblim==10.0);%index of blood vessel
                    if (~isempty(bldidx))
                        lblim(bldidx)=4.0;
                    end
                    inflidx = find(core_lblim==11.0);
                    if (~isempty(inflidx))
                        lblim(inflidx)=5.0;
                    end
                    nervidx = find(core_lblim==12.0);
                    if (~isempty(nervidx))
                        lblim(nervidx)=6.0;
                    end
                    corpidx = find(core_lblim==8.0);
                    if (~isempty(corpidx))
                        lblim(corpidx)=7.0;
                    end
                end
                writeTIFF(lblim,strcat(curpath,label_file_name),'int8');
              
                
            else
                lblim=imread(strcat(curpath,label_file_name));
            end
            figure(1);
            imagesc(lblim);
        end
    end
    close(h);

end