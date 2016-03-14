 clc;
 clear all;
 close all;
ntextons = 51;
%Generate the data and train a classifier
[filenames,glandnames]=findFileName('F:\QPI data\');
data_dir='F:\QPI data\texdir\';
label_dir='F:\QPI data\label\';
svm_dir = 'F:\QPI data\svm\ws20\';
total_sample_num=4*67000;
sample_per_im = 2000;
nfilestotal = 0;
h = waitbar(0,'Reading data textons...');
data = zeros(0,ntextons);
label=zeros(0,1);
filelist=cell(0,1);
testmode = 2; 
%Testing mode: 1: test on the validation subsample data
%2: test on the full-size image
%3: compute the Performance curve of the SVM algorithm
if (~exist('svmstruct_ws80.mat','file'))
    if (~exist('sampleddata_ws80.mat','file'))
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
                    tx_hist_file_name = strcat(label_name,'_texton_hist_80.mat');
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
                    curlabel=[ones(sample_per_im,1);(-1)*ones(sample_per_im,1)]; %Gland 1, stroma -1
                    curdata=[gland_hist1;stroma_hist1];
                    data(end+1:end+2*sample_per_im,:)=curdata;
                    label(end+1:end+2*sample_per_im)=curlabel;
                    filelist{end+1,1}=label_name;

               end
        end
        save('sampleddata_ws80.mat','data','label','filelist');
    else
        load('sampleddata_ws80.mat');
    end
    disp('Begin trainng process...')
    testset = data(1:2*end/3,:);
    testlabel = label(1:2*end/3,:);
    %Downsample the data and do SVM training
    opts = statset('Display','iter','MaxIter',8000000,'TolX',0.01);
    svmtruct=svmtrain(data,label,'method','SMO','options',opts,'Kernel_Function','rbf');
    save 'svmstruct_ws80.mat' svmtruct;

else
    disp('Loading trained classifier...')
    load('svmstruct_ws40.mat');
    disp('Loading validation data...');
    switch testmode
        case 1
            load('sampleddata_ws40.mat','data','label');
            valset = data(round(2*end/3)+1:end,:);
            vallabel = label(:,round(2*end/3)+1:end);
            nsampleperbatch = 1000;
            nbatch=floor(size(valset,1)/nsampleperbatch);
            nsampletest = nsampleperbatch*nbatch;
            valset = valset(1:nsampletest,:);
            vallabel = vallabel(1,1:nsampletest);
            outputval=zeros(nsampletest,1);
            for batchidx=1:nbatch
                disp(['Batch idx ', num2str(batchidx)]);
                outputval((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:)=...
                    svmclassify(svmtruct,valset((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:));
                outputlabel = outputval((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:);
                groundtruth = vallabel(1,(batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch);
                accu_rate = (length(find(outputlabel(:)==groundtruth(:))))/length(groundtruth)
                disp(['Cross-validation accuracy' num2str(accu_rate)]);
            end
            %Downsample the data and do SVM training
            accu_rate = (length(find(outputval(:)==vallabel(:)))/length(vallabel));
            disp(['Cross-validation accuracy' num2str(accu_rate)]);
        case 2
            
            for classidx=1:4 %Go through different classes
               nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
               for sampleIdx=1:nsamples
                    nfilestotal =nfilestotal+1; 
                    nfilestotal =nfilestotal+1; 
                    waitbar(sampleIdx/nsamples,h,'Progress...')
                    cur_file_name = filenames{classidx,1}{sampleIdx,1};
                    %Check to see if the label file exist
                    dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                    slash_pos = strfind(cur_file_name,'\');
                    label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                    disp(['Adding data of: ' label_name ' ...']);
                    tx_hist_file_name = strcat(label_name,'_texton_hist_80.mat');
                    svm_file_name = strcat(label_name,'_svm_200k_ws80.mat');
                    if (~exist(strcat(svm_dir,svm_file_name)))
                        load(strcat(data_dir,tx_hist_file_name),'histim');
                        load(strcat(label_dir,label_name,'_resized.mat'));

                        nrows = 512;
                        ncols = 512;
                        %Validation with downsampled image
                        newhistim = zeros(nrows,ncols,ntextons);
                        lblim=imresize(lblim,[nrows ncols],'nearest');
                        fim = (-3)*ones(nrows,ncols);
                        %Resize the labeling so that the stroma goes to -1
                        tempidx= find(lblim==2);
                        lblim(tempidx)=-1; %Readjust the label of the stroma

                        %Load the value of a classifier to compute
                        shiftval = svmtruct.ScaleData.shift;
                        scaleval = svmtruct.ScaleData.scaleFactor;
                        supvect = svmtruct.SupportVectors;%Get the support vectors
                        alpha = svmtruct.Alpha; %Note that alpha is positive for the 1st group and -1 for the second group
                        bias = svmtruct.Bias;
                        kerfunc = svmtruct.KernelFunction;
                        kerfuncargs = svmtruct.KernelFunctionArgs;

                        testdata = zeros(nrows*ncols,ntextons);
                        for bandidx=1:ntextons
                            bandidx;
                            newhistim(:,:,bandidx)=imresize(histim(:,:,bandidx),[nrows,ncols],'nearest');
                        end
                        histim = newhistim;

                        %Line up to convert from 3D data to 2D data
                        for bandidx=1:ntextons
                             curband = histim(:,:,bandidx);
                             curvect = curband(:);
                             testdata(:,bandidx)=curvect(:);
                        end
                        clear histim;
                        clear newhistim;
                        %Now, just take classify gland vs stroma
                        tempidx = find(lblim~=0);
                        vallabel = lblim(tempidx);
                        valset = testdata(tempidx,:);

                        nsampleperbatch =4096;
                        nbatch=ceil(size(valset,1)/nsampleperbatch);
                        outputval=zeros(size(vallabel));
                        f=zeros(size(vallabel)); %This is the value of the classification function

                        for batchidx=1:nbatch
                            disp(['Batch index: ', num2str(batchidx), '/' num2str(nbatch)]);
                            curevalset =  valset((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,size(valset,1)),:);                       
                            curevalset = bsxfun(@plus,curevalset,shiftval);%Data normalization
                            curevalset = bsxfun(@times,curevalset,scaleval);
                            temp_out=-(kerfunc(supvect,curevalset,kerfuncargs{:})'*alpha(:) + bias(:));
                            f((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,size(valset,1)),:)=...
                                temp_out;
                            temp_label = (temp_out>=0)+(-1)*(temp_out<0);
                            %temp_label=svmclassify(svmtruct,curevalset);
                            outputval((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,size(valset,1)),:) = temp_label;
                            outputlabel = temp_label;%Current output label
                            groundtruth = vallabel((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,size(valset,1)));
                            accu_rate = (length(find(outputlabel(:)==groundtruth(:))))/length(groundtruth)
                            disp(['Patch accuracy: ' num2str(accu_rate)]);

                        end
                        accu_rate = (length(find(outputval(:)==vallabel(:)))/length(vallabel));
                        disp(['Cross-validation accuracy: ' num2str(accu_rate)]);

                        %Display an image after classification
                        tempidx = find(lblim~=0);
                        newlblim = lblim;
                        newlblim(tempidx)=outputval;
                        groundtruth = lblim(tempidx);
                        fim(tempidx)=f;
                        figure(1);
                        subplot(121);
                        imagesc(reshape(newlblim,[nrows ncols]));
                        title('SVM results');
                        subplot(122);
                        imagesc(reshape(lblim,[nrows ncols]));
                        title('Ground truth');
                        figure(2);
                        imagesc(fim);
                        colorbar;
                        title('f');
                        [xg,yg]=perfcurve(groundtruth,f,1);
                        figure(3);
                        plot(xg,yg)
                        xlabel('False positive rate (gland)'); ylabel('True positive rate (gland)');
                        [xs,ys]=perfcurve(groundtruth,-f,-1);
                        figure(4);
                        plot(xs,ys,'r');
                        xlabel('False positive rate (stroma)'); ylabel('True positive rate (stroma)');

                       % save(strcat(svm_dir,svm_file_name),'fim','newlblim','f','outputval','lblim',...);
                       %     'xs','xg','ys','yg');
                       save(svm_file_name,'fim','newlblim','f','outputval','lblim',...);
                            'xs','xg','ys','yg','tempidx');
                   end
               end
            end
      case 3 %ROC curve on the training set...
            
            %Compute the performance curve by using the function f
            %First, compute the slack variables which will measure the
            %degree of mis-classifications...
            load('sampleddata.mat','data','label');
            trainset = data(1:2*end/3,:);
            
            shiftval = svmtruct.ScaleData.shift;
            scaleval = svmtruct.ScaleData.scaleFactor;
            trainset = bsxfun(@plus,trainset,shiftval);%Data normalization
            trainset = bsxfun(@times,trainset,scaleval);
            
            supvect = svmtruct.SupportVectors;%Get the support vectors
            alpha = svmtruct.Alpha; %Note that alpha is positive for the 1st group and -1 for the second group
            bias = svmtruct.Bias;
            kerfunc = svmtruct.KernelFunction;
            kerfuncargs = svmtruct.KernelFunctionArgs;
            
            nsampleperbatch = 500;
            nbatch=floor(size(trainset,1)/nsampleperbatch);
            nsampletest = nsampleperbatch*nbatch;
            trainset = trainset(1:nsampletest,:);
            trainlabel = label(1,1:nsampletest);
            f=zeros(nsampletest,1);
            
            for batchidx=1:nbatch
                disp(['Batch idx ', num2str(batchidx)]);
                f((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:)=...
                    -kerfunc(supvect,trainset((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:),kerfuncargs{:})'*alpha(:) + bias(:);
            end
            
            %f is the score for the +1 class
            [x,y]=perfcurve(trainlabel(1:2500),f(1:2500),1);
            plot(x,y);
            xlabel('False positive rate'); ylabel('True positive rate');
    end
end


