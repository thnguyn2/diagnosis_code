function roc_texton    
clc;
    clear all;
    close all;
    labelpath = 'H:\QPI data\label\';
    datapath = 'H:\QPI data\';
    svm_dir1 = 'H:\QPI data\svm\texton_ws20_512\';
    svm_dir2 = 'H:\QPI data\svm\texton_ws40_512\';
    svm_dir3 = 'H:\QPI data\svm\texton_ws80_512\';
    [filenames,glandnames]=findFileName(datapath);
    
    %part1: compute the distribution for each image and save it
    h = waitbar(0,'Calculating oriented energy progress...');
    nsamples_per_data = 4000;
    nradius = 3;
    npixels = 512^2;
    nrows = 512;
    ncols = 512;
    
    score_arr_svm = zeros(0,nradius);
    label_arr_svm = zeros(0,nradius);
    newlabelvec = zeros(0,1);
    gtlabelvec = zeros(0,1);
    pixlabelidx=1;
    
    xs_svm = cell(nradius,1);
    ys_svm = cell(nradius,1);
    auc_svm = cell(nradius,1);
    sampleidx=1;
    colorarr='rbcykmg'; 
    
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       accu =zeros(2,1);
       for sampleIdx=1:nsamples
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            cur_file_name = 'H:\QPI data\33D1.tif';
            label_file_name = strcat(label_name,'.mat');
            org_im = imread(cur_file_name); %Load the original data
            load(strcat(labelpath,label_file_name));
            org_im = imresize(org_im,[2048 2048],'nearest');
            lblim = imresize(lblim,[2048 2048],'nearest');
            org_im = cast(org_im,'single');
            figure(5);
            imagesc(org_im/65536*4.0-0.5);
            colorbar;
            colormap jet;
            axis off;
            
            load(strcat(labelpath,label_file_name));
            lblim = imresize(lblim,[nrows ncols],'nearest');
            glandidx = find(lblim==1);
            stromaidx = find(lblim==2);
            %Compute ROC curve for other radius
            pixelidx=[1:npixels];
            [y_coord,x_coord]=ind2sub([nrows ncols], pixelidx); %Row and column indices of current pixels
            y_coord = y_coord(:);
            x_coord = x_coord(:);
            
            
            svm_file_name1 = strcat(label_name,'_svm_200k_ws20.mat'); 
            load(strcat(svm_dir1,svm_file_name1));
            disp(['Calculating phase histogram: ' label_name ' ...']);                     
            %Add the output energy for evaluation
            glandsvm1 = fim(glandidx); %Compute the score for gland and stroma
            stromasvm1 = fim(stromaidx);
            glandsvm1 = glandsvm1(1:2*floor(size(glandsvm1,1)/nsamples_per_data):end);
            stromasvm1 = stromasvm1(1:2*floor(size(stromasvm1,1)/nsamples_per_data):end);
            score = [glandsvm1(:);stromasvm1(:)];
            labeldata=[ones(size(glandsvm1));2*ones(size(stromasvm1))];
            label_arr_svm(sampleidx:sampleidx + length(score)-1,1)=labeldata;
            score_arr_svm(sampleidx:sampleidx + length(score)-1,1)=score;                 
            [xs_svm{1},ys_svm{1},t,auc_svm{1}]=perfcurve(label_arr_svm(:,1),score_arr_svm(:,1),1);
            
            svm_file_name2 = strcat(label_name,'_svm_200k_ws40.mat'); 
            load(strcat(svm_dir2,svm_file_name2));
            %Add the output energy for evaluation
            glandsvm2 = fim(glandidx); %Compute the score for gland and stroma
            stromasvm2 = fim(stromaidx);
            glandsvm2 = glandsvm2(1:2*floor(size(glandsvm2,1)/nsamples_per_data):end);
            stromasvm2 = stromasvm2(1:2*floor(size(stromasvm2,1)/nsamples_per_data):end);
            score = [glandsvm2(:);stromasvm2(:)];
            labeldata=[ones(size(glandsvm2));2*ones(size(stromasvm2))];
            label_arr_svm(sampleidx:sampleidx + length(score)-1,2)=labeldata;
            score_arr_svm(sampleidx:sampleidx + length(score)-1,2)=score;                 
            [xs_svm{2},ys_svm{2},t,auc_svm{2}]=perfcurve(label_arr_svm(:,2),score_arr_svm(:,2),1);
            
            %Obtain a segmentaion map for different images with the best
            %texton feature
            newlblim = zeros(size(lblim));
            thresh = -0.0295;
            candidateidx = find(lblim~=0);
            candidatelabel = (-1)*ones(size(candidateidx));
            glandidx = find(fim(candidateidx)>thresh);
            candidatelabel(glandidx)=1;
            newlblim(candidateidx)=candidatelabel;
            lumenidx = find(newlblim==0);
            newlblim(lumenidx)=-1;
            %Apply image dilation
            se = strel('disk',1);
            idxim = im2bw(newlblim,0);
            idxim = imclose(idxim,se);
            newlblim = cast(idxim,'single')*2.0-1.0;
            newlblim(lumenidx)=0;
            colorbar
            
            [newlblim] = procseg(newlblim); %Another processing step to remove sall are
            figure(3);
            imagesc(newlblim);
            title('Labeling results from texton');
            axis off;
            figure(4);
            imagesc(lblim);
            title('Human labeling');
            drawnow;
            axis off;
            
            disp('Current confusion matrix...');
            m = computefusionmat(lblim(candidateidx), newlblim(candidateidx));
           
            accu(1)=(accu(1)*(sampleIdx-1)+m(1,1))/sampleIdx
            accu(2)=(accu(2)*(sampleIdx-1)+m(2,2))/sampleIdx
            
            
            newlabelvec(pixlabelidx:pixlabelidx + length(lblim(candidateidx))-1,1) =...
                newlblim(candidateidx);
            gtlabelvec(pixlabelidx:pixlabelidx + length(lblim(candidateidx))-1,1) =...
                lblim(candidateidx);
            pixlabelidx = pixlabelidx + length(lblim(candidateidx));
            %gtlabelvec = zeros(0,1);
           
           
            
            
            
            
            svm_file_name3 = strcat(label_name,'_svm_200k_ws80.mat'); 
            load(strcat(svm_dir3,svm_file_name3));
            %Add the output energy for evaluation
            glandsvm3 = fim(glandidx); %Compute the score for gland and stroma
            stromasvm3 = fim(stromaidx);
            glandsvm3 = glandsvm3(1:2*floor(size(glandsvm3,1)/nsamples_per_data):end);
            stromasvm3 = stromasvm3(1:2*floor(size(stromasvm3,1)/nsamples_per_data):end);
            score = [glandsvm3(:);stromasvm3(:)];
            labeldata=[ones(size(glandsvm3));2*ones(size(stromasvm3))];
            label_arr_svm(sampleidx:sampleidx + length(score)-1,3)=labeldata;
            score_arr_svm(sampleidx:sampleidx + length(score)-1,3)=score;                 
            [xs_svm{3},ys_svm{3},t,auc_svm{3}]=perfcurve(label_arr_svm(:,3),score_arr_svm(:,3),1);
            
            
            figure(1);
            plot(xs_svm{1},ys_svm{1},colorarr(1),'linewidth',2);
            hold on;
            plot(xs_svm{2},ys_svm{2},colorarr(2),'linewidth',2);
            hold on;
            plot(xs_svm{3},ys_svm{3},colorarr(3),'linewidth',2);
            
            drawnow;       
            sampleidx = sampleidx + length(score);
            figure(1);
            hold off;
            legend('20','40','80');
            title('ROC curve for the mean phase');
            auc_svm{2}
            
            
         
  
          end

 %   save('textonroc.mat','xs_svm','ys_svm','auc_svm','label_arr_svm','score_arr_svm');
    close(h);
    end
end

function m = computefusionmat(gtlabel,newlabel)
    m=zeros(2,2);
    npixels = length(gtlabel(:));
    glandidx = find(gtlabel==1);
    nglands = length(glandidx(:));
    stromaidx = find(gtlabel==-1);
    nstromas = length(stromaidx(:));
    glandmatch = intersect(find(gtlabel==1),find(newlabel==1));
    m(1,1)=length(glandmatch)/nglands;
    m(1,2)=1-m(1,1);
    stromamatch = intersect(find(gtlabel==-1),find(newlabel==-1));
    m(2,2)=length(stromamatch)/nstromas;
    m(2,1)=1-m(2,2);
    disp([num2str(m(1,1)) num2str(m(1,2))]);
    disp([num2str(m(2,1)) num2str(m(2,2))]);
    
end
