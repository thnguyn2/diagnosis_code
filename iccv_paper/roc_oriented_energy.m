     clc;
    clear all;
    close all;
    labelpath = 'H:\QPI data\label\';
    datapath = 'H:\QPI data\';
    histpath = 'H:\QPI data\text_dir_hist\';
    oepath = 'H:\QPI data\texdir\'
    [filenames,glandnames]=findFileName(datapath);
    
    %part1: compute the distribution for each image and save it
    h = waitbar(0,'Calculating oriented energy progress...');
    nsamples_per_data = 4000;
    nradius = 5;
    npixels = 2048^2;
    nrows = 2048;
    ncols = 2048;
    score_arr_oe = zeros(0,nradius);
    label_arr_oe = zeros(0,nradius);
    score_arr_hod = zeros(0,1);
    label_arr_hod = zeros(0,1);
    xs_oe = cell(nradius,1);
    ys_oe = cell(nradius,1);
    auc_oe = cell(nradius,1);
    sampleidx=1;
     colorarr='rbcykmg';  
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            label_file_name = strcat(label_name,'.mat');
            fr_file_name = strcat(label_name,'_lm_fr.mat');
            texture_dir_hist_file_name =  strcat(label_name,'_text_dir_hist.mat');
            load(strcat(histpath,texture_dir_hist_file_name),'klsym');

            load(strcat(oepath,fr_file_name),'fe','fo','emag','emap','edir');
            
            disp(['Calculating phase histogram: ' label_name ' ...']);
            org_im = imread(cur_file_name); %Load the original data
            load(strcat(labelpath,label_file_name));
            org_im = imresize(org_im,[2048 2048],'nearest');
            lblim = imresize(lblim,[2048 2048],'nearest');
            org_im = cast(org_im,'single');
            glandidxorg = find(lblim==1);
            stromaidxorg = find(lblim==2);
            
            %Compute ROC curve for other radius
            pixelidx=[1:npixels];
            [y_coord,x_coord]=ind2sub([2048 2048], pixelidx); %Row and column indices of current pixels
            y_coord = y_coord(:);
            x_coord = x_coord(:);
            
            oe = fe+fo; %Compute the oriented energy
           
            for radidx = 1:nradius

                %Add the output energy for evaluation
                curoe = oe(:,:,radidx);
                glandoe = curoe(glandidxorg); %Compute the score for gland and stroma
                stromaoe = curoe(stromaidxorg);
                glandoe = glandoe(1:2*floor(size(glandoe,1)/nsamples_per_data):end);
                stromaoe = stromaoe(1:2*floor(size(stromaoe,1)/nsamples_per_data):end);

                score = [glandoe(:);stromaoe(:)];
                labeldata=[ones(size(glandoe));2*ones(size(stromaoe))];
                label_arr_oe(sampleidx:sampleidx + length(score)-1,radidx)=labeldata;
                score_arr_oe(sampleidx:sampleidx + length(score)-1,radidx)=score;
               
                
                
                [xs_oe{radidx},ys_oe{radidx},t,auc_oe{radidx}]=perfcurve(label_arr_oe(:,radidx),score_arr_oe(:,radidx),2);
                figure(1);
                plot(xs_oe{radidx},ys_oe{radidx},colorarr(radidx),'linewidth',2);
                hold on;
                
            end
            
            %Compute the feature for kl distance
            glandhod = klsym(glandidxorg); %Compute the score for gland and stroma
            stromahod = klsym(stromaidxorg);
            glandhod = glandhod(1:2*floor(size(glandhod,1)/nsamples_per_data):end);
            stromahod = stromahod(1:2*floor(size(stromahod,1)/nsamples_per_data):end);

            scorehod = [glandhod(:);stromahod(:)];
            labeldata=[ones(size(glandhod));2*ones(size(stromahod))];
            label_arr_hod(sampleidx:sampleidx + length(score)-1,1)=labeldata;
            score_arr_hod(sampleidx:sampleidx + length(score)-1,1)=scorehod;
               
                
                
            [xs_hod,ys_hod,t,auc_hod]=perfcurve(label_arr_hod(:,1),score_arr_hod(:,1),2);
            plot(xs_hod,ys_hod,colorarr(1),'linewidth',2);
            hold on;

            
            
            sampleidx = sampleidx + length(score);
            figure(1);
            hold off;
            legend('0.50 um','0.71 um','1.01 um','1.43 um','2.02 um','hod');
            title('ROC curve for the mean phase');
          end
    end
 
    save('oehod.mat','xs_oe','ys_oe','auc_oe','label_arr_oe','score_arr_oe','label_arr_hod','score_arr_hod',...
        'xs_hod','ys_hod','auc_hod');

    close(h);
