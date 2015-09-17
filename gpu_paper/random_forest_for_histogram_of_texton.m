 clc;
 clear all;
 close all;
ntextons = 51;
%Generate the data and train a classifier
addpath(strcat(cd(cd('../..')),'\RF_code'));
if (~exist('forest_Nov_1.mat'))
    [filenames,glandnames]=findFileName('H:\QPI data\');
    data_dir='H:\QPI data\texdir\';
    label_dir='H:\QPI data\label\';
    total_sample_num=4*67000;
    sample_per_im = 2000;
    nfilestotal = 0;
    h = waitbar(0,'Reading data textons...');
    data = zeros(0,ntextons);
    label=zeros(0,1);
    filelist=cell(0,1);
    if (~exist('sampleddata.mat','file'))
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
                    disp(['Adding data of: ' label_name ' ...']);
                    tx_hist_file_name = strcat(label_name,'_texton_hist_40.mat');
                    load(strcat(data_dir,tx_hist_file_name),'histim');
                    load(strcat(label_dir,label_name,'_resized.mat'));
                    glandidx = find(lblim==1);
                    stromaidx = find(lblim==2);
                    lumenidx = find(lblim==0);
                    ngland = length(glandidx);
                    nstroma = length(stromaidx);
                    gland_hist = zeros(ngland,ntextons);
                    stroma_hist = zeros(nstroma,ntextons);

                    for bandidx=1:ntextons
                        curband = histim(:,:,bandidx);
                        glandvect=curband(glandidx);
                        stromavect=curband(stromaidx);
                        lumenvect=curband(lumenidx);
                        gland_hist(:,bandidx)=glandvect(:);
                        stroma_hist(:,bandidx)=stromavect(:);
                    end
                    glandidx=randperm(ngland);
                    stromaidx=randperm(nstroma);
                    gland_hist1=gland_hist(glandidx(1:sample_per_im),:);
                    stroma_hist1=stroma_hist(stromaidx(1:sample_per_im),:);
                    curlabel=[ones(sample_per_im,1);(-1)*ones(sample_per_im,1)];
                    curdata=[gland_hist1;stroma_hist1];
                    data(end+1:end+2*sample_per_im,:)=curdata;
                    label(end+1:end+2*sample_per_im)=curlabel;
                    filelist{end+1,1}=label_name;

               end
        end
        save('sampleddata.mat','data','label','filelist');
    else
        load('sampleddata.mat');
    end
    label=label(:);
    traindata = data(1:round(2*end/3),:);
    trainlabel = label(1:round(2*end/3),:);

    testdata=data(round(2*end/3)+1:end,:);
    testlabel = label(round(2*end/3)+1:end,:);
    params.ntrees=1;
    params.max_level=20;
    params.min_node_size=5; %Minimum node size
    params.eval_oob=1;      %Evaluate the Out of bag error
    params.var_impt=1;      %Evaluate the importance of variable
    nfeatures=size(traindata,2);
    ntrain=size(traindata,1);
    [forest,label_arr,oob_idx]=class_RF_train(traindata,trainlabel,params);
    save('forest_Nov_1.mat','forest','label_arr','oob_idx');
else
    load('forest_Nov_1.mat');
    load('sampleddata.mat');
    testdata=data(round(2*end/3)+1:end,:);
    testlabel = label(round(2*end/3)+1:end,:);
    
    %[outputestlabel,p_hat]=eval_RF(testdata,forest,label_arr);
    [outputestlabel,p_hat]=eval_RF(testdata(89321,:),forest,label_arr);
    
    idx=find(outputestlabel~=testlabel);
   disp(['Prediction error: ' num2str(length(idx)/length(testlabel))]);

end

% clear histim
% load('H:\QPI data\texdir\33D1_texton_hist_40.mat','histim');
% load('H:\QPI data\label\33D1_resized.mat');
% %Validation with downsampled image
% newhistim = zeros(256,256,ntextons);
% lblim=imresize(lblim,[256 256],'nearest');
% 
% testdata = zeros(256^2,ntextons);
% for bandidx=1:ntextons
%     bandidx;
%     newhistim(:,:,bandidx)=imresize(histim(:,:,bandidx),[256,256]);
% end
% histim = newhistim;
% a=histim(28,88,:);
% 
% for bandidx=1:ntextons
%     bandidx
%     curband = histim(:,:,bandidx);
%     curvect = curband(:);
%     testdata(:,bandidx)=curvect(:);
% end
% testlabel=lblim(:);
% [outputestlabel,p_hat]=eval_RF(testdata,forest,label_arr);   
% %Display the probability of class
% chan1 = reshape(p_hat(1,:),[256 256]);
% figure(1)
% imagesc(chan1);
% idx=find(outputestlabel~=testlabel);
% disp(['Prediction error: ' num2str(length(idx)/length(testlabel))]);
% 
