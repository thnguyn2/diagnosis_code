function diagnosis_demo
    clc;
    clear all;
    close all;
    do_cancer_vs_noncancer=0;
    do_cancer_only =1;
    
    addpath(strcat(cd(cd('..')),'\support'));
    datapath = 'G:\TMA_cores_and_diagnosis\';
    diagpath = 'G:\TMA_cores_and_diagnosis\diag\';
    disp('Finding the list of all labels');
    [filenames,glandnames,gradelist]=findFileNameFromROIs(datapath);
    save ('filename.mat','filenames','glandnames','gradelist');
    mode = 'gt'; %Two different kinds of modes 'gt' and 'real' for using two different kind of data (ground truth and the real data)
    %Train or update morphological features
    stromabandrad = 45;
    updateGlandSize(filenames,strcat(datapath,'label\'),strcat(datapath,'diag\'),strcat(datapath,'texdir\'),mode,stromabandrad);%Compute the diagnosis information and save it
    
    
    load(strcat(diagpath,strcat('glandmorp_',mode,'_',num2str(stromabandrad),'.mat')),'feat','dataclass');
    %Note that the first 25 features are for gland's morphological features
    %the next 120 is the histogram of stroma
    %the final 120 is the histogram of glands
    
    %Find the first and the last sample index for each class
    ngrades = size(filenames,1);
    sampleidx = zeros(ngrades,2);
    nsamples_added=0;
    %Create a matrix containing the first and the last sample indices for
    %each class
    for classidx=1:ngrades %Go through different classes
       sampleidx(classidx,1)=nsamples_added+1;
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            label_file_name = strcat(label_name,'_resized.mat');
            if (exist(label_file_name,'file'))               
                nsamples_added = nsamples_added + 1;
            end   
       end
       sampleidx(classidx,2)=nsamples_added;
    end

    nsamplesperclass=sampleidx(:,2)-sampleidx(:,1)+1;
    ntextons=50;
    
    class23idx=find(strcmp(gradelist,'23')==1);
    class32idx=find(strcmp(gradelist,'32')==1);
    class33idx=find(strcmp(gradelist,'33')==1);
    class34idx=find(strcmp(gradelist,'34')==1);
    class35idx=find(strcmp(gradelist,'35')==1);
    class43idx=find(strcmp(gradelist,'43')==1);
    class44idx=find(strcmp(gradelist,'44')==1);
    class45idx=find(strcmp(gradelist,'45')==1);
    class53idx=find(strcmp(gradelist,'53')==1);
    class54idx=find(strcmp(gradelist,'54')==1);
    class55idx=find(strcmp(gradelist,'55')==1);
    classBPidx=find(strcmp(gradelist,'BP')==1);
    classHGidx=find(strcmp(gradelist,'HG')==1);
    classNMidx=find(strcmp(gradelist,'NM')==1);
    
    idx23 = sampleidx(class23idx,1):sampleidx(class23idx,2);idx23=idx23(:);
    idx32 = sampleidx(class32idx,1):sampleidx(class32idx,2);idx32=idx32(:);    
    idx33 = sampleidx(class33idx,1):sampleidx(class33idx,2);idx33=idx33(:);
    idx34 = sampleidx(class34idx,1):sampleidx(class34idx,2);idx34=idx34(:);    
    idx35 = sampleidx(class35idx,1):sampleidx(class35idx,2);idx35=idx35(:);    
    idx43 = sampleidx(class43idx,1):sampleidx(class43idx,2);idx43=idx43(:);
    idx44 = sampleidx(class44idx,1):sampleidx(class44idx,2);idx44=idx44(:);
    idx45 = sampleidx(class45idx,1):sampleidx(class45idx,2);idx45=idx45(:);
    idx53 = sampleidx(class53idx,1):sampleidx(class53idx,2);idx53=idx53(:);
    idx54 = sampleidx(class54idx,1):sampleidx(class54idx,2);idx54=idx54(:);
    idx55 = sampleidx(class55idx,1):sampleidx(class55idx,2);idx55=idx55(:);
  
    idxHG = sampleidx(classHGidx,1):sampleidx(classHGidx,2);idxHG=idxHG(:);
    idxBP = sampleidx(classBPidx,1):sampleidx(classBPidx,2);idxBP=idxBP(:);
    idxNM = sampleidx(classNMidx,1):sampleidx(classNMidx,2);idxNM=idxNM(:);
    
     %Part 2: Test with multiple class classifier
     data23 = feat(idx23,:);
     data32 = feat(idx32,:);
     data33 = feat(idx33,:);
     data34 = feat(idx34,:);
     data35 = feat(idx35,:);
     data43 = feat(idx43,:);
     data44 = feat(idx44,:);
     data45 = feat(idx45,:);
     data53 = feat(idx53,:);
     data54 = feat(idx54,:);
     data55 = feat(idx55,:);       
     
     dataHG = feat(idxHG,:);
     dataBP = feat(idxBP,:);
     dataNM = feat(idxNM,:);
     data =[data23;data32;data33;...
            data34;data43;...
            data35;data53;data44;data45;data54;data55;...
            dataHG;dataBP;dataNM];
     %For cancer/non-cancer classification, just go with the most important
     %features. RF tells that the feature 25th is the most important one
     fulldata=data;
     data = data(:,26:end); %Get the histogram of stroma only
     
     species = dataclass([idx23;idx32;idx33;...
                         idx34;idx43;...
                         idx35;idx53;idx44;idx45;idx54;idx55;...
                         idxHG;idxBP;idxNM],1);
    if (do_cancer_vs_noncancer)
         cannormspecies = [ones(size([idx23;idx32;idx33;...
                             idx34;idx43;...
                             idx35;idx53;idx44;idx45;idx54;idx55]));...
                            zeros(size([idxHG;idxBP;idxNM]))];%Cancer/normal class

         canceridx = find(cannormspecies==1);             
         normidx=find(cannormspecies==0);
         ncancer = length(canceridx);
         nnormal = length(normidx);
         figure(1);
         subplot(211);
         %Draw the histogram of stroma
         for i=1:ncancer
             plot(data(canceridx(i),:),'r');
             hold on;
             drawnow;
         end
         title('Cancer case');
         subplot(212);
         for j=1:nnormal
             plot(data(normidx(j),:),'b');
             hold on;
             drawnow;
         end
         title('Non-cancer')

         %Find and draw 2 most important features
         [idx, z] = rankfeatures(data',cannormspecies,'Criterion','roc');
         [zshort,V]=sort(z,'descend');             
         figure(2);
         gscatter(data(:,V(1)),data(:,V(2)),cannormspecies,'rgbcy');   
         xlabel(num2str(V(1)));
         ylabel(num2str(V(2)));
         legend('Cancer','Benign');

         nfeatures = 10;
         data =  data(:,V(1:nfeatures));
         %Verify with the location of basal cells

         disp('[Part 1]: Cancer /non-cancer classifier');
         kval = 10;
         nround = 50; %Number of round
         cftotal = zeros(2,2);
         for roundidx=1:nround %Each round start with a different random seeds
             rng(roundidx*5)
             CVO = cvpartition(cannormspecies,'kfold',kval); %Go with 4 since we have at least 1 sample for each class
             err = zeros(CVO.NumTestSets,1);
             cp_arr=zeros(CVO.NumTestSets,1);
             cfmat = cell(CVO.NumTestSets,1);
             classifiertype = 2;
             accum_varimp = zeros(1,size(data,2));
             for i = 1:CVO.NumTestSets
                  trIdx = CVO.training(i);
                  teIdx = CVO.test(i);


                  switch classifiertype
                      case 1
                          %Bayesian classifier
                          c = NaiveBayes.fit(data(trIdx,:),cannormspecies(trIdx,:));
                          yout = c.predict(data(teIdx,:));
                      case 2
                          %Svm classifier - works for data after PCA reduction
                          %only since the number of observation is smaller than
                          %the number of features in our problem
                          svmstruct = svmtrain(data(trIdx,:),cannormspecies(trIdx,:),'kernel_function','rbf');
                          yout = svmclassify(svmstruct,data(teIdx,:));

                      case 3
                          ntrees = 100;
                          rfstruct = RTrees;
                          rfstruct.train(data(trIdx,:),cannormspecies(trIdx,:),'MaxNumOfTreesInTheForest',ntrees,...
                            'ForestAccuracy',0.05,'MaxDepth',8,'Prior',[10 1],'VarType','Categorical','CalcVarImportance',1); %Putting more weights on stroma.
                          yout = rfstruct.predict(data(teIdx,:));
                          evalfeatimp=0;
                          if (evalfeatimp)
                              %Evaluate the importance of features
                              varimp = rfstruct.getVarImportance();
                              accum_varimp = (accum_varimp*(i-1)+varimp)/i; %Accumulated feature importance
                              [U,V]=sort(accum_varimp,'descend');
                              feat1idx=V(1);feat2idx=V(2);
                              gscatter(data(:,feat1idx),data(:,feat2idx),cannormspecies,'rgbcy');
                              xlabel(num2str(feat1idx));
                              ylabel(num2str(feat2idx));
                              drawnow
                          end
                      case 4 %Discriminant analysis. Assume each class has a Gaussian distribution with its own mean and covariance - It will give the model for each class
                          classstruct =  ClassificationDiscriminant.fit(data(trIdx,:),cannormspecies(trIdx,:),'discrimType','diagQuadratic');
                          [yout score cost] = predict(classstruct,data(teIdx,:));

                      case 5 %Nearest neighbot classifier
                                mdl = ClassificationKNN.fit(data(trIdx,:),cannormspecies(trIdx,:),'NumNeighbors',3);
                                %L = loss(mdl,cancerdata(teIdx,:),cancerspecies(teIdx,:));
                                %disp(['K-NN lost: ' num2str(L)]);
                                yout = predict(mdl,data(teIdx,:)); %Predict the label of the output


                  end
                  cp = classperf(cannormspecies(teIdx,:),yout);
                  cp_arr(i)=cp.CorrectRate;
                  confusionmat(cannormspecies(teIdx,:),yout);
                  cfmat{i,1} = confusionmat(cannormspecies(teIdx,:),yout);
                  err(i) = sum(yout~=cannormspecies(teIdx))/length(yout);
             end
            curcf = zeros(size(cfmat{1,1}));
            for cfmatidx=1:size(cfmat,1)
                 curcf = curcf + cfmat{cfmatidx,1};
            end
            curcf = curcf/kval;
            cftotal = (cftotal*(roundidx-1)+curcf)/roundidx;
            cftotalnew =cftotal./repmat(sum(cftotal,2),[1, size(cftotal,2)])*100

            %disp(['Round: ' num2str(roundidx)]);
         end 
         disp(['Error after ' num2str(nround) ' with' num2str(kval) '-fold cross-validation']);
         cftotalnew =cftotal./repmat(sum(cftotal,2),[1, size(cftotal,2)])*100;
 
    end
     
     
     %========================================================================
     disp('Step 2: classification for the non-cancer class- BPH & HGPIN');
     data = [dataHG;dataBP;dataNM];
    % data = [dataHG;dataBP];
     data = data(:,1:end); %Get the histogram of stroma only
     normspecies = [zeros(size([idxHG]));ones(size([idxBP]));ones(size(idxNM))];
    % normspecies = [zeros(size([idxHG]));ones(size([idxBP]))];
     hgidx = find(normspecies==0);             
     bpidx=find(normspecies==1);
     nmidx =find(normspecies==2);
%      nhg = length(hgidx);
%      nbp = length(bpidx);
%      nnm = length(nmidx);
%      [idx, z] = rankfeatures(data',normspecies,'Criterion','roc');
%      [zshort,V]=sort(z,'descend');             
%      figure(2);
%      gscatter(data(:,V(1)),data(:,V(2)),normspecies,'rgbcy');   
%      xlabel(num2str(V(1)));
%      ylabel(num2str(V(2)));
%      legend('BPH','HGPIN');
     
%     nfeatures = 2;
%     data =  data(:,V(1:nfeatures));
     %Verify with the location of basal cells
      
     disp('[Part 2]: BPH /HGPIN/Normal classifier');
     kval = 9;
     nround = 50; %Number of round
     ncategories = length(unique(normspecies));
     cftotal = zeros(ncategories,ncategories);
    
     for roundidx=1:nround %Each round start with a different random seeds
         %rng(roundidx*20)
         CVO = cvpartition(normspecies,'kfold',kval); %Go with 4 since we have at least 1 sample for each class
         err = zeros(CVO.NumTestSets,1);
         cp_arr=zeros(CVO.NumTestSets,1);
         cfmat = cell(CVO.NumTestSets,1);
         classifiertype = 4;
         accum_varimp = zeros(1,size(data,2));
         for i = 1:CVO.NumTestSets
              trIdx = CVO.training(i);
              teIdx = CVO.test(i);
              
                         
              switch classifiertype
                  case 1
                      %Bayesian classifier
                      c = NaiveBayes.fit(data(trIdx,:),normspecies(trIdx,:));
                      yout = c.predict(data(teIdx,:));
                  case 2
                      %Svm classifier - works for data after PCA reduction
                      %only since the number of observation is smaller than
                      %the number of features in our problem
                      svmstruct = svmtrain(data(trIdx,:),normspecies(trIdx,:),'kernel_function','rbf');
                      yout = svmclassify(svmstruct,data(teIdx,:));
                      
                  case 3
                      ntrees = 100;
                      rfstruct = RTrees;
                      rfstruct.train(data(trIdx,:),normspecies(trIdx,:),'MaxNumOfTreesInTheForest',ntrees,...
                        'ForestAccuracy',0.05,'MaxDepth',8,'Prior',[1 1],'VarType','Categorical','CalcVarImportance',1); %Putting more weights on stroma.
                      yout = rfstruct.predict(data(teIdx,:));
                      evalfeatimp=0;
                      if (evalfeatimp)
                          %Evaluate the importance of features
                          varimp = rfstruct.getVarImportance();
                          accum_varimp = (accum_varimp*(i-1)+varimp)/i; %Accumulated feature importance
                          [U,V]=sort(accum_varimp,'descend');
                          feat1idx=V(1);feat2idx=V(2);
                          gscatter(data(:,feat1idx),data(:,feat2idx),cannormspecies,'rgbcy');
                          xlabel(num2str(feat1idx));
                          ylabel(num2str(feat2idx));
                          drawnow
                      end
                  case 4 %Discriminant analysis. Assume each class has a Gaussian distribution with its own mean and covariance - It will give the model for each class
                      classstruct =  ClassificationDiscriminant.fit(data(trIdx,:),normspecies(trIdx,:),'discrimType','diagQuadratic');
                      [yout score cost] = predict(classstruct,data(teIdx,:));
                            
                  case 5 %Nearest neighbot classifier
                            mdl = ClassificationKNN.fit(data(trIdx,:),normspecies(trIdx,:),'NumNeighbors',3);
                            %L = loss(mdl,cancerdata(teIdx,:),cancerspecies(teIdx,:));
                            %disp(['K-NN lost: ' num2str(L)]);
                            yout = predict(mdl,data(teIdx,:)); %Predict the label of the output
      
                      
              end
              cp = classperf(normspecies(teIdx,:),yout);
              cp_arr(i)=cp.CorrectRate;
              confusionmat(normspecies(teIdx,:),yout);
              cfmat{i,1} = confusionmat(normspecies(teIdx,:),yout);
              err(i) = sum(yout~=normspecies(teIdx))/length(yout);
         end
        curcf = zeros(size(cfmat{1,1}));
        for cfmatidx=1:size(cfmat,1)
             curcf = curcf + cfmat{cfmatidx,1};
        end
        curcf = curcf/kval;
        cftotal = (cftotal*(roundidx-1)+curcf)/roundidx;
        %disp(['Round: ' num2str(roundidx)]);
        cftotalnew =cftotal./repmat(sum(cftotal,2),[1, size(cftotal,2)])*100
     end 
     disp(['Error after ' num2str(nround) ' with' num2str(kval) '-fold cross-validation']);
     cftotalnew

  
       
     
   
     
     %====================================================================
     disp('Step 3: Classification for the cancer class');
     datalt7 = [data23;data32;data33];
     dataeq7 = [data34;data43];
     datagt7 = [data35;data53;data44;data45;data54;data55];
     cancerdata = [datalt7;dataeq7;datagt7];

     idxlt7 = [idx23;idx32;idx33];
     idxeq7 = [idx34;idx43];
     idxgt7 = [idx35;idx53;idx44;idx45;idx54;idx55];

     cancerspecies = [6*ones(size(idxlt7));7*ones(size([idxeq7]));7*ones(size(idxgt7))];
     
%     [idx, z] = rankfeatures(cancerdata',cancerspecies,'Criterion','roc');
%     [zshort,V]=sort(z,'descend');
%       figure(2);
%      gscatter(cancerdata(:,V(1)),cancerdata(:,V(2)),cancerspecies,'rgbcy');   
%      xlabel(num2str(V(1)));
%      ylabel(num2str(V(2)));
%  
%     nfeatures = 2;
%     selectedfeat = V(1:nfeatures);
%     cancerdata=cancerdata(:,selectedfeat);
     
     kval = 10;
     nround = 50; %Number of round
     ncategories = length(unique(cancerspecies));
     cftotal = zeros(ncategories,ncategories);
     %1fold cross-validation may perform worse than training error
     for roundidx=1:nround %Each round start with a different random seeds  
             CVO = cvpartition(cancerspecies,'kfold',kval); %Go with 4 since we have at least 1 sample for each class
             err = zeros(CVO.NumTestSets,1);
             cp_arr=zeros(CVO.NumTestSets,1);
             cfmat = cell(CVO.NumTestSets,1);
             classifiertype = 1;
             accum_varimp = zeros(1,size(data,2));
             for i = 1:CVO.NumTestSets
                 trIdx = CVO.training(i);
                 teIdx = CVO.test(i);
                    switch classifiertype
                     case 1
                          %Bayesian classifier
                          c = NaiveBayes.fit(cancerdata(trIdx,:),cancerspecies(trIdx,:));
                          yout = c.predict(cancerdata(teIdx,:));
                     case 2
                          %Svm classifier - works for data after PCA reduction
                          %only since the number of observation is smaller than
                          %the number of features in our problem
                          svmstruct = svmtrain(cancerdata(trIdx,:),cancerspecies(trIdx,:),'kernel_function','rbf');
                          yout = svmclassify(svmstruct,cancerdata(teIdx,:));

                     case 3
                          ntrees = 100;
                          rfstruct = RTrees;
                          rfstruct.train(cancerdata(trIdx,:),cancerspecies(trIdx,:),'MaxNumOfTreesInTheForest',ntrees,...
                               'ForestAccuracy',0.05,'MaxDepth',8,'Prior',[1 1],'VarType','Categorical','CalcVarImportance',1); %Putting more weights on stroma.
                          yout = rfstruct.predict(cancerdata(teIdx,:));
                          evalfeatimp=0;
                          if (evalfeatimp)
                              %Evaluate the importance of features
                              varimp = rfstruct.getVarImportance();
                              accum_varimp = (accum_varimp*(i-1)+varimp)/i; %Accumulated feature importance
                              [U,V]=sort(accum_varimp,'descend');
                              feat1idx=V(1);feat2idx=V(2);
                              gscatter(cancerdata(:,feat1idx),cancerdata(:,feat2idx),cancerspecies,'rgbcy');
                              xlabel(num2str(feat1idx));
                              ylabel(num2str(feat2idx));
                              drawnow
                          end               
                      case 4 %Discriminant analysis
                         %classstruct =  ClassificationDiscriminant.fit(cancerdata(trIdx,:),cancerspecies(trIdx,:),'discrimType','diagQuadratic');
                         %[yout score cost] = predict(classstruct,cancerdata(teIdx,:));
                         yout = classify(cancerdata(teIdx,:),cancerdata(trIdx,:),cancerspecies(trIdx,:),'quadratic');
                         
                      case 5 %Nearest neighbot classifier
                            mdl = ClassificationKNN.fit(cancerdata(trIdx,:),cancerspecies(trIdx,:),'NumNeighbors',3);
                            %L = loss(mdl,cancerdata(teIdx,:),cancerspecies(teIdx,:));
                            %disp(['K-NN lost: ' num2str(L)]);
                            yout = predict(mdl,cancerdata(teIdx,:)); %Predict the label of the output
                     end
                      cp = classperf(cancerspecies(teIdx,:),yout);
                      cp_arr(i)=cp.CorrectRate;
                      cfmat{i,1} = confusionmat(cancerspecies(teIdx,:),yout);
                      err(i) = sum(yout~=cancerspecies(teIdx))/length(yout);
             end
             curcf = zeros(size(cfmat{1,1}));
             for cfmatidx=1:size(cfmat,1)
                     curcf = curcf + cfmat{cfmatidx,1};
             end
             curcf = curcf/kval;
             cftotal = (cftotal*(roundidx-1)+curcf)/roundidx;
             cftotalnew =cftotal./repmat(sum(cftotal,2),[1, size(cftotal,2)])*100;
     end
     disp(['[Part 3]: Error after ' num2str(nround) ' with' num2str(kval) '-fold cross-validation']);
     cftotalnew   
     
     
    %Generate the decision boundary if tge data is 2D    
    if (size(cancerdata,2)==2) %Draw the decision domain
        [x,y]=meshgrid(linspace(min(cancerdata(:,1)),max(cancerdata(:,1)),100),linspace(min(cancerdata(:,2)),max(cancerdata(:,2))));
        x=x(:);
        y=y(:);
         switch classifiertype
               case 1
                          %Bayesian classifier
                          c = NaiveBayes.fit(cancerdata,cancerspecies);
                          yout = c.predict(cancerdata);
                          disp('Mis-classification for the training');
                          confusionmat(cancerspecies,yout)
                          yout = c.predict([x y]);
               case 2
                          svmstruct = svmtrain(cancerdata,cancerspecies,'kernel_function','rbf');
                          %Compute the missclassification error
                          yout = svmclassify(svmstruct,cancerdata);
                          disp('Mis-classification for the training');
                          confusionmat(cancerspecies,yout)
                          yout = svmclassify(svmstruct,[x y]);
                          
                case 3
                          ntrees = 100;
                          rfstruct = RTrees;
                          rfstruct.train(cancerdata,cancerspecies,'MaxNumOfTreesInTheForest',ntrees,...
                               'ForestAccuracy',0.05,'MaxDepth',8,'Prior',[1 1],'VarType','Categorical','CalcVarImportance',1); %Putting more weights on stroma.
                          yout = rfstruct.predict(cancerdata);
                          disp('Mis-classification for the training');
                          confusionmat(cancerspecies,yout)
                          yout = rfstruct.predict([x y]);                                                   
                      
               case 4
                          yout = classify(cancerdata,cancerdata,cancerspecies,'quadratic');
                          disp('Mis-classification for the training');
                          confusionmat(cancerspecies,yout)
                          yout = classify([x,y],cancerdata,cancerspecies,'quadratic');
                   
               case 5 %Nearest neighbot classifier
                            mdl = ClassificationKNN.fit(cancerdata,cancerspecies,'NumNeighbors',1);
                            yout = predict(mdl,cancerdata); %Predict the label of the output
                            
                            confusionmat(cancerspecies,yout)
                            yout = predict(mdl,[x y]); %Predict the label of the output
                            
         end                      
         figure(3);
         gscatter(x,y,yout,'rbgcy','.',2);
         hold on;
         gscatter(cancerdata(:,1),cancerdata(:,2),cancerspecies,'rbgcy');
         xlabel(strcat('Feat ',num2str(V(1))));
         ylabel(strcat('Feat ',num2str(V(2))));
    
                    
    else
    end

end

function [mcsvmstruct,uniclass]=multisvmtrain(X,Y)
    %This function defines the multi-class SVM. Each row of X is an
    %observation
    uniclass= unique(Y);
    nclass=length(uniclass);
    mcsvmstruct = cell(nclass,1);
    for classidx=1:nclass
        %disp(['Training for class: ' num2str(classidx)]);
        class = ((Y==uniclass(classidx))-0.5)*2;
        opts = statset('Display','off','MaxIter',20000,'TolX',0.01);
        mcsvmstruct{classidx,1}=svmtrain(X,class,'method','SMO','options',opts,'Kernel_Function','linear');        
    end   
end

function [outclass,rank]=multisvmclassify(multisvm,valX,uniclass)
    %Perform SVM classification and assign the class to X based on the
    %classifier that gives max functional distance
    nclassf=size(multisvm,1); %Number of SVM classifier
    rank=zeros(size(valX,1),nclassf);%Output of the classifier
    for classidx=1:nclassf
          %Load the value of a classifier to compute
          shiftval = multisvm{classidx,1}.ScaleData.shift;
          scaleval = multisvm{classidx,1}.ScaleData.scaleFactor;
          supvect = multisvm{classidx,1}.SupportVectors;%Get the support vectors
          alpha = multisvm{classidx,1}.Alpha; %Note that alpha is positive for the 1st group and -1 for the second group
          bias = multisvm{classidx,1}.Bias;
          kerfunc = multisvm{classidx,1}.KernelFunction;
          kerfuncargs = multisvm{classidx,1}.KernelFunctionArgs;
          f=zeros(size(valX,1),1);                              
          curevalset = bsxfun(@plus,valX,shiftval);%Data normalization
          curevalset = bsxfun(@times,curevalset,scaleval);
          temp_out=-(kerfunc(supvect,curevalset,kerfuncargs{:})'*alpha(:) + bias(:));
          rank(:,classidx)=temp_out;                 
    end   
    [maxval,outclassidx]=max(rank,[],2);
    outclass=uniclass(outclassidx);
    
end

function error=classerror(xtrain,ytrain,xtest,ytest)
    yout = classify(xtest,ytrain,ytrain);
    %Try NaiveBayes classifier
%     c = NaiveBayes.fit(xtrain,ytrain);
%     yout = c.predict(xtest);
     error = sum(~(yout~=ytest));%Use the mistmatch as the gold standard for evalutating the feature. 
end


function classnum=cell2number(curcell)
    %Convert the cell to classnumber '+1'-> 1, '-1'>-1
    nsamples=size(curcell,1);
    classnum = zeros(size(curcell));
    for sampleidx=1:nsamples
        classnum(sampleidx,1)=str2double(cell2mat(curcell(sampleidx,1)));
    end
end
