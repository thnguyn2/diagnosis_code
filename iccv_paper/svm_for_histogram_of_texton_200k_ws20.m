
import cv.*
clc;
 clear all;
 close all;
%Train a classifier with for cores with the glands
%Last update" Feb 28th 2015
ntextons = 50;
classf_type = 1;
%Generate the data and train a classifier
datapath = 'E:\Dino_data\TMA_cores_and_diagnosis\';
addpath(strcat(cd(cd('..')),'\support'));
[filenames,glandnames,classname,tiffiles]=findFileNameFromROIs(datapath);
svm_dir = strcat(datapath,'diag\'); %Path to store the data for SVM classifier

radius = 60;
nrows = 3072;
ncols = 3072;


testmode = 2; 
ngrades = size(filenames,1);
computingplatform = 1;
retrain=0;
reeval = 1;
nclasses = 2;
if ((~exist(strcat(svm_dir,strcat('rfstruct_ws',num2str(radius),'_March_18th_2015_',num2str(nclasses),'_classes.mat')),'file'))||(retrain==1))
    recreatesampledata=0;
    if ((~exist(strcat(svm_dir,strcat('sampleddata_March_9th_2015_ws',num2str(radius),'.mat')),'file'))||...
            (recreatesampledata==1))
        total_sample_num=0;
        sample_per_im = 8000;%This is the number of gland pixels and the number of stroma pixels for each image with groundtruth
        nfilestotal = 0;

        data = zeros(0,ntextons);
        label=zeros(0,1);
        filelist=cell(0,1);
        %h = waitbar(0,'Reading data textons...');
        texton_dir=strcat(datapath,'texdir\');
        label_dir=strcat(datapath,'label\');
        for classidx=1:ngrades %Go through different classes
               nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
               for sampleIdx=1:nsamples
                    nfilestotal =nfilestotal+1; 
                   % waitbar(sampleIdx/nsamples,h,'Progress...')
                    cur_file_name = filenames{classidx,1}{sampleIdx,1};
                    %Check to see if the label file exist
                    dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                    slash_pos = strfind(cur_file_name,'\');
                    label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                    disp(['Adding data of: ' label_name ' ...']);
                    tx_hist_file_name = strcat(label_name,'_texton_hist_',num2str(radius),'.mat');
                    load(strcat(texton_dir,tx_hist_file_name));
                    %Load the label image
                    lbl_file_name = strcat(label_name,'_label.tif');
                    lblim=imread(strcat(label_dir,lbl_file_name));
                    lblim = imresize(lblim,[size(histim,1) size(histim,2)],'nearest');
                    figure(1);
                    imagesc(lblim);
                    drawnow;
                    %Make sure we have all the pixels for the thin layer of
                    %stroma
                    stromaim = (lblim==2);
                    stromaimerode = imerode(stromaim,strel('disk',30));
                    stromaim = stromaim-stromaimerode;
                    stromaidxmoreimp = find(stromaim==1);
                   
                    stromaidx=stromaidxmoreimp;
                    
                    glandidx = find(lblim==1);
                    lumenidx = find(lblim==0);
                    scrtidx = find(lblim==3);
                    bldidx = find(lblim==4);
                    inflidx = find(lblim==5);
                    nervidx = find(lblim==6);
                    corpidx = find(lblim==7);
                    
                    nlumen = length(lumenidx);
                    ngland = length(glandidx);
                    nstroma = length(stromaidx);
                    nscrt=length(scrtidx);
                    ninfl = length(inflidx);
                    nbld = length(bldidx);
                    nnerv = length(nervidx);
                    ncorp = length(corpidx);
                    
                    gland_hist = zeros(ngland,ntextons);
                    stroma_hist = zeros(nstroma,ntextons);
                    lumen_hist = zeros(nlumen,ntextons);
                    if (nscrt>0)
                        scrt_hist = zeros(nscrt,ntextons);
                    end
                    if (ninfl>0)
                        infl_hist = zeros(ninfl,ntextons);
                    end
                    if (nbld>0)
                        bld_hist =zeros(nbld,ntextons);
                    end
                    if (nnerv>0)
                        nerv_hist = zeros(nnerv,ntextons);
                    end
                    if (ncorp>0)
                        corp_hist = zeros(ncorp,ntextons);
                    end
                    for bandidx=1:ntextons
                        curband = histim(:,:,bandidx);
                        glandvect=curband(glandidx);
                        stromavect=curband(stromaidx);
                        lumenvect=curband(lumenidx);
                        if (nscrt>0)
                            scrtvect=curband(scrtidx);
                            scrt_hist(:,bandidx)=scrtvect(:);
                        end
                        if (ninfl>0)
                            inflvect=curband(inflidx);
                            infl_hist(:,bandidx)=inflvect(:);
                        end
                        if (nbld>0)
                            bldvect = curband(bldidx);
                            bld_hist(:,bandidx)=bldvect(:);
                        end
                        if (nnerv>0)
                            nervect = curband(nervidx);
                            nerv_hist(:,bandidx)=nervect(:);
                        end
                        if (ncorp>0)
                            corpvect = curband(corpidx);
                            corp_hist(:,bandidx)=corpvect(:);
                        end
                        
                        gland_hist(:,bandidx)=glandvect(:);
                        stroma_hist(:,bandidx)=stromavect(:);
                        lumen_hist(:,bandidx)=lumenvect(:);           
                                       
                       
                        
                        
                    end
                    glandidx=randperm(ngland);
                    stromaidx=randperm(nstroma);
                    lumenidx = randperm(nlumen);
                    if (nscrt>0)
                        scrtidx = randperm(nscrt);
                    end
                    if (ninfl>0)
                        inflidx = randperm(ninfl);
                    end
                    if (nbld>0)
                        bldidx = randperm(nbld);
                    end
                    if (nnerv>0)
                        nervidx = randperm(nnerv);
                    end
                    if (ncorp>0)
                        corpidx = randperm(ncorp);
                    end
                    
                    if (ngland>sample_per_im)
                        ngland = sample_per_im;
                    end
                    if (ngland>0)                        
                        gland_hist1=gland_hist(glandidx(1:ngland),:);
                    end
                    if (nstroma>sample_per_im)
                        nstroma = sample_per_im; 
                    end
                    if (nstroma>0)
                        stroma_hist1=stroma_hist(stromaidx(1:nstroma),:);
                    end
                    if (nlumen>sample_per_im)
                        nlumen = sample_per_im;
                    end
                    if (nlumen>0)                        
                        lumen_hist1=lumen_hist(lumenidx(1:nlumen),:);
                    end
                   
                    minorclass_data=zeros(0,ntextons);
                    minorclass_data(end+1:end+ngland,:)=gland_hist1;
                    minorclass_data(end+1:end+nstroma,:)=stroma_hist1;
                    minorclass_data(end+1:end+nlumen,:)=lumen_hist1;
                    
                    if (nscrt>sample_per_im)
                        nscrt = sample_per_im;
                    end
                    if (nscrt>0)
                        scrt_hist1 = scrt_hist(scrtidx(1:nscrt),:);                 
                        minorclass_data(end+1:end+nscrt,:)=scrt_hist1;
                    end
                    if (ninfl>sample_per_im)
                        ninfl = sample_per_im;
                    end
                    if (ninfl>0)
                        infl_hist1 = infl_hist(inflidx(1:ninfl),:);
                        minorclass_data(end+1:end+ninfl,:)=infl_hist1;
                    end
                    if (nbld>sample_per_im)
                        nbld = sample_per_im;
                    end
                    if (nbld>0)
                        bld_hist1 = bld_hist(bldidx(1:nbld),:);
                        minorclass_data(end+1:end+nbld,:)=bld_hist1;
                    end
                    if (ncorp>sample_per_im)
                        ncorp = sample_per_im;
                    end
                    if (ncorp>0)
                        corp_hist1 = corp_hist(corpidx(1:ncorp),:);
                        minorclass_data(end+1:end+ncorp,:)=corp_hist1;
                    end
                    if (nnerv>sample_per_im)
                        nnerv = sample_per_im;
                    end
                    if (nnerv>0)
                        nerv_hist1 = nerv_hist(nervidx(1:nnerv),:);
                        minorclass_data(end+1:end+nnerv,:)=nerv_hist1;
                    end
                    curlabel=[ones(ngland,1);2*ones(nstroma,1);0*ones(nlumen,1);...
                        3*ones(nscrt,1);4*ones(nbld,1);5*ones(ninfl,1);6*ones(nnerv,1);7*ones(ncorp,1)]; %Gland 1, stroma 2; Lumen 0; scretion 3; blood 4; inflammatory cell 5
                    %nerv: 6, corpora amylacel:7
                  
                    curdata= minorclass_data;   
                    data(end+1:end+length(curlabel(:)),:)=curdata;
                    label(end+1:end+length(curlabel(:)))=curlabel;
                    filelist{end+1,1}=label_name;
                    total_sample_num = total_sample_num + length(curlabel);
               end
        %       close(h);
        end
        data = cast(data,'single');
        save(strcat(svm_dir,'sampleddata_March_9th_2015_ws',num2str(radius),'.mat'),'data','label','filelist','-v7.3');
    else
        load(strcat(svm_dir,'sampleddata_March_9th_2015_ws',num2str(radius),'.mat'));
    end
    disp('Begin trainng process...')
    %Just keep a few classes of interest, eliminate others
    keepidx = find(label>=1);
    data = data(keepidx,:);
    label=label(keepidx);
    outlieridx = find(label>=3);
    label(outlieridx)=2; %The outlier are also consider as a part of stroma
    nsamples = size(data,1);
    ntrees = 50;
    fullidx=1:nsamples;
    trainidx=find(mod(fullidx,2)~=0);
    testidx=find(mod(fullidx,2)==0);
    data = data(1:nsamples,:);
    label = label(1:nsamples);
    traindata = data(trainidx,:);
    trainlabel = label(trainidx);
    
   
    testdata = data(testidx,:);
    testlabel = label(testidx);
   % save(strcat(svm_dir,'sampleddata_Feb_28th_2015_ws',num2str(radius),'.mat'),'traindata','trainlabel','testdata','testlabel');
    tic
    if ((computingplatform==0)) %MATLAB
        if (classf_type==0) %SVM classifier
            %Downsample the data and do SVM training
            opts = statset('Display','iter','MaxIter',8000000,'TolX',0.01);
            svmtruct=svmtrain(traindata,trainlabel,'method','SMO','options',opts,'Kernel_Function','rbf');
            save(strcat(svm_dir,'svmstruct_ws',num2str(radius),'.mat'),svmtruct);
        else
            matlabpool open; %Use parallel computing toolbox for faster computation
            %Random Forest classifiers
            opt = statset('UseParallel','always');
            rfstruct = TreeBagger(ntrees,traindata,trainlabel,'OOBPred','on','Options',opt,...
                'Method','classification','NVarToSample',11,'NPrint',1,'MinLeaf',5);
            plot(oobError(rfstruct));
            matlabpool close;
            save(strcat(svm_dir,'rfstruct_ws',num2str(radius),'_Feb_28_2015_3_classes.mat'),'rfstruct');
            

        end
    else 
        if (classf_type==0) %SVM
            %Fast but not reliable
            %Not good enough for our data since the number of support vector
            %must be large. In OpenCV, no more than 2000 support vectors are
            %allows
            svmstruct = SVM;
            trainingres = svmstruct.train(traindata,trainlabel,'KernelType','RBF','C',0.1);
            svmstruct.save(strcat(svm_dir,'svmstruct_ws',num2str(radius),'.mat'))
        else %RF
            rfstruct = RTrees;
            rfstruct.train(traindata,trainlabel,'MaxNumOfTreesInTheForest',ntrees,...
            'ForestAccuracy',0.05,'MaxDepth',25,'Prior',[1 20],'VarType','Categorical','NActiveVars',11); %Putting more weights on stroma.
            rfstruct.save(strcat(svm_dir,'rfstruct_ws',num2str(radius),'_March_18th_2015_',num2str(nclasses),'_classes.mat'))
       end
    end
    timeElap = toc;
    disp(['Training time: ' num2str(timeElap/60.0), ' (min)']);
    
else
    %Evaluation phase----
    %Testing mode:
    %1: test on the validation subsample data
    %2: test on the full-size image
    %3: compute the Performance curve of the SVM algorithm
    reeval=1;
    if (classf_type==0) %SVM classifier
        if (computingplatform==0) %MATLAB
            load(strcat(svm_dir,'svmstruct_ws',num2str(radius),'.mat'))
        else
        end
    else %Random Forest classifer
        if (computingplatform==0) %MATLAB
            load(strcat(svm_dir,'rfstruct_ws',num2str(radius),'.mat'))
        else      
            if (reeval)
                rfstruct = RTrees;
                rfstruct.load(strcat(svm_dir,'rfstruct_ws',num2str(radius),'_March_18th_2015_',num2str(nclasses),'_classes.mat'));
            end

        end
    end
    disp('Loading trained classifier...')
    if (reeval)
       % load(strcat(svm_dir,'train&testdata_Aug_12th_ws',num2str(radius),'.mat'));
    end
    disp('Loading validation data...');
    reeval=1;
    switch testmode
        case 1
            valset = testdata;
            vallabel = testlabel;
            if (computingplatform==0)
                if (classf_type==0) %SVM classifier
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
                else %RF classifier
                    nsampleperbatch = 1000;
                    nbatch=floor(size(valset,1)/nsampleperbatch);
                    nsampletest = nsampleperbatch*nbatch;
                    valset = valset(1:nsampletest,:);
                    vallabel = vallabel(1,1:nsampletest);
                    outputval=zeros(nsampletest,1);
                    for batchidx=1:nbatch
                        %disp(['Batch idx ', num2str(batchidx)]);
                        cures=predict(rfstruct,valset((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:));
                        outputval((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:)=cell2number(cures);
                        outputlabel = outputval((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:);
                        groundtruth = vallabel(1,(batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch);
                        accu_rate = (length(find(outputlabel(:)==groundtruth(:))))/length(groundtruth);
                        %disp(['Cross-validation accuracy: ' num2str(accu_rate)]);
                    end
                end
            else %OpenCVMEX computation
                if (classf_type==1) %RF classifier
                    nsampleperbatch = 80000;
                    nbatch=floor(size(valset,1)/nsampleperbatch);
                    nsampletest = nsampleperbatch*nbatch;
                    valset = valset(1:nsampletest,:);
                    vallabel = vallabel(1,1:nsampletest);
                    outputval=zeros(nsampletest,1);
                    for batchidx=1:nbatch
                        disp(['Batch idx ', num2str(batchidx)]);
                        cures=rfstruct.predict(valset((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:));
                        outputval((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:)=cures;
                        outputlabel = outputval((batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch,:);
                        groundtruth = vallabel(1,(batchidx-1)*nsampleperbatch+1:batchidx*nsampleperbatch);
                        accu_rate = (length(find(outputlabel(:)==groundtruth(:))))/length(groundtruth);
                        %disp(['Cross-validation accuracy: ' num2str(accu_rate)]);
                    end
                end

            end
            %Downsample the data and do SVM training
            accu_rate = (length(find(outputval(:)==vallabel(:)))/length(vallabel));
            disp(['Cross-validation accuracy: ' num2str(accu_rate)]);
            c=confusionmat(vallabel(:),outputval)
      case 2
            data_dir = strcat(datapath,'texdir\');
            nfilestotal =0;
            ngrades = size(filenames,1);
            nsamples = size(tiffiles,1);
           
            for sampleIdx=1:nsamples
           % for classidx=1:ngrades %Go through different classes
           %    nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
           %    for sampleIdx=1:nsamples
                    nfilestotal =nfilestotal+1; 
                    %waitbar(sampleIdx/nsamples,h,'Progress...')
                    cur_file_name = tiffiles{sampleIdx,1};
                    %cur_file_name = filenames{classidx,1}{sampleIdx,1};
                    %Check to see if the label file exist
                    dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
                    slash_pos = strfind(cur_file_name,'\');
                    label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
                    disp(['Adding data of: ' label_name ' ...']);
                    tx_hist_file_name = strcat(label_name,'_texton_hist_',num2str(radius),'.mat');
                    if (computingplatform==0) %This area may have been depricated
                        svm_file_name = strcat(label_name,'_svm_200k_ws',num2str(radius),'.mat');
                        if (~exist(svm_file_name))
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

                            nsampleperbatch =1000;
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
                    else
                         %Random Forest with OpenCVMex
                         fuzzy_output = 1;
                        
                         if (fuzzy_output==0)   
                            seg_im_name = strcat(cur_file_name(1:end-4),'_seg_ws',num2str(radius),'.tif');
                         else
                            seg_im_name = strcat('H:\TMA_cores_and_diagnosis\label\',num2str(nclasses),'class_results\',label_name,'_',num2str(radius),...
                                 '_',num2str(nclasses),'cl.tif');  
                            seg_im_name_bf = strcat(cur_file_name(1:end-4),'_seg_ws',num2str(radius),'_fuzzy_bf.tif');  
                         end
                         if ((~exist(seg_im_name,'file'))||(reeval==1))
                             load(strcat(data_dir,tx_hist_file_name));
                             nrows = size(histim,1);
                             ncols = size(histim,2);
                             ntextons = size(histim,3);
                             testdata = zeros(nrows*ncols,ntextons);
                             %Line up to convert from 3D data to 2D data
                             for bandidx=1:ntextons
                                 curband = histim(:,:,bandidx);
                                 curvect = curband(:);
                                 testdata(:,bandidx)=curvect(:);
                             end
                             clear histim;
                             %if we don't have the label image, generate it
                             %now
                             gtavail = 0;
                             org_im =imread(cur_file_name);
                             lblim = findLabelImg(org_im);
                             
                             %Now, just take classify gland vs stroma
                             tempidx = find(lblim~=0);
                             if (gtavail)
                                vallabel = lblim(tempidx);                                
                             end
                             valset = testdata(tempidx,:);
                             nevalsample = size(valset,1);
                             nsampleperbatch =80000;
                             nbatch=ceil(size(valset,1)/nsampleperbatch);
                             outputlabel=zeros(size(tempidx));
                             outputprob=zeros(size(tempidx));
                           
                             
                             for batchidx=1:nbatch
                                disp(['Batch index: ', num2str(batchidx), '/' num2str(nbatch)]);                              
                                if (fuzzy_output==0)
                                    cures=rfstruct.predict(valset((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,nevalsample),:));
                                    outputval((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,nevalsample),:) = cures;
                                    if (gtavail) %If we have groundtruth, then compute the accuracy
                                        groundtruth = vallabel((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,nevalsample));
                                        accu_rate = (length(find(cures(:)==groundtruth(:))))/length(groundtruth)
                                        disp(['Patch accuracy: ' num2str(accu_rate)]);
                                    end
                                else
                                    if (nclasses>=3)
                                    %This is with more than 3 classes
                                        cures=rfstruct.predict(valset((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,nevalsample),:));
                                        outputlabel((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,nevalsample),:)=cures;
                                    else    
                                        %Two classes, produce a probability
                                        %map though
                                        fuzzy_cures=rfstruct.predict_prob(valset((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,nevalsample),:));                                
                                        outputprob((batchidx-1)*nsampleperbatch+1:min(batchidx*nsampleperbatch,nevalsample),:) = fuzzy_cures;                                 
                               
                                    end
                                end
                                
                             end
                             if ((fuzzy_output==0)&(gtavail))
                                accu_rate = (length(find(outputval(:)==vallabel(:)))/length(vallabel));
                                disp(['Cross-validation accuracy: ' num2str(accu_rate)]);
                             end

                             %Display an image after classification
                             tempidx = find(lblim~=0);
                             newlblim = lblim;
                             if (nclasses>=3)
                                newlblim(tempidx)=outputlabel;
                             else
                                newlblim(tempidx)=outputprob;
                             end
                             groundtruth = lblim(tempidx);
                             figure(1);
                             imagesc(reshape(newlblim,[nrows ncols]));
                             title('RF results');
                             drawnow;
                             if (nclasses>=3)
                                writeTIFF(cast(newlblim,'uint8'),strcat('H:\TMA_cores_and_diagnosis\label\',num2str(nclasses),'class_results\',label_name,'_',num2str(radius),...
                                 '_',num2str(nclasses),'cl.tif'),'int8');
                             else
                                 writeTIFF(newlblim,strcat('H:\TMA_cores_and_diagnosis\label\',num2str(nclasses),'class_results\',label_name,'_',num2str(radius),...
                                 '_',num2str(nclasses),'cl.tif'));
                            
                             end
         
                         else
                         end
                         if ((fuzzy_output)||(~exist(seg_im_name_bf))) %if we go with the fuzzy output and the bilateral filter
                             %seg_im_name = imread(seg_im_name);
                             
                         end
            
                    end
              % end
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


