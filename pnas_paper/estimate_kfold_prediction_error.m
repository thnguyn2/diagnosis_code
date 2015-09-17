function [accu,dev,mdl,varimportance] = estimate_kfold_prediction_error(data,label,k,type)
    %Esmate the prediction error based on k-fold cross validation
    %Inputs:
    %   data, label: a matrix for the input data and label data for
    %   classification. Each row is an observation...
    %   k: number of data cluster that we used for cross-validation
    %   type: type of the classifer
    %Outputs:
    %   accu, dev: mean and standard deviation of the estimation...
    
    label_arr = unique(label);
    nclass = length(label_arr(:));
    %Get the number of samples per class
    splperclass = zeros(nclass,1);
    indexperclass = cell(nclass,1);
    varimportance = 0;
    for classidx=1:nclass
        splperclass(classidx,1)=sum(label==label_arr(classidx));
        temp = find(label==label_arr(classidx));
        indexperclass{classidx,1}=temp;
    end
    accu = 0;
    dev = 0;
    accu_mat = zeros(k,1);
    if (k>min(splperclass))
        return;
    else
        varimportance = zeros(1,size(data,2));
        for kidx=1:k
                 valfeat = zeros(0,size(data,2));
                 vallabel = zeros(0,1);
                 validxused = zeros(0,1);
                 trainfeat = zeros(0,size(data,2));
                 trainlabel = zeros(0,1);
                 %Get the k-portion to validation
                 for classidx=1:nclass
                    curnsplinclass=length(indexperclass{classidx,1}(:)); 
                    validx = indexperclass{classidx,1}((kidx-1)*floor(curnsplinclass/k)+1:...
                       kidx*floor(curnsplinclass/k));
                    validxused(end+1:end+length(validx),:)=validx;
                    valfeat(end+1:end+length(validx),:)=data(validx,:);%Add to the validation set    
                    vallabel(end+1:end+length(validx),:)=label(validx,:);
                 end
                 %Get the unused index for training...
                 trainidx = setdiff(1:size(data,1),validxused);
                 trainfeat(end+1:end+length(trainidx),:)=data(trainidx,:);
                 trainlabel(end+1:end+length(trainidx),:)=label(trainidx,:);
                 %Now do the classification...
                 switch type
                     case {'knn'}
                          k_neigh =1; %Number of nearest neighbors
                          mdl = ClassificationKNN.fit(trainfeat,trainlabel,'NumNeighbors',k_neigh);
                          predlabel = predict(mdl,valfeat);
                          %disp(['Confusion matrix at fold ' num2str(kidx) '-th']);
                          %cm = confusionmat(vallabel,predlabel)%Compute the confusion matrix...
                          %Compute the accuracy mat...
                          accu_mat(kidx,1)=mean(vallabel==predlabel);
                          disp(['Accuracy: ' num2str(accu_mat(kidx,1))]);
                     case {'svm'}
                          if (nclass>2)
                              error('SVM does not work when the number of class is larger than 2');
                              return;
                          else
                              opts = statset('Display','iter','MaxIter',16000000,'TolX',0.01);
                              mdl=svmtrain(trainfeat,trainlabel,'method','SMO','options',opts,'Kernel_Function','rbf');
                              predlabel=svmclassify(mdl,valfeat,'showplot','true');
                              disp(['Confusion matrix at fold ' num2str(kidx) '-th']);
                              cm = confusionmat(vallabel,predlabel)%Compute the confusion matrix...
                              %Compute the accuracy mat...
                              accu_mat(kidx,1)=mean(vallabel==predlabel);
                              disp(['Accuracy: ' num2str(accu_mat(kidx,1))]);
                          end
                    case {'classtree'}
                        mdl = classregtree(trainfeat,num2str(trainlabel));
                        predlabel = eval(mdl,valfeat);
                        disp(['Confusion matrix at fold ' num2str(kidx) '-th']);
                        cm = confusionmat(vallabel,str2num(cell2mat(predlabel)))%Compute the confusion matrix...
                        accu_mat(kidx,1)=mean(vallabel==str2num(cell2mat(predlabel)));
                        disp(['Accuracy: ' num2str(accu_mat(kidx,1))]);
                    %Compute a tree baggin classifier
                     case {'treebagger'}
                         ntrees = 200;
                         mdl = TreeBagger(ntrees,trainfeat,trainlabel,'OOBPred','on','OOBVarImp','on');
                          %Plot the out-of-bagg error...
                         %plot(oobError(mdl));
                         %xlabel('Number of trees...');
                         %ylabel('OOB classification errors...');
                         %Predict the output label from the tree bagger...
                         predlabel = predict(mdl,valfeat);
                         disp(['Confusion matrix at fold ' num2str(kidx) '-th']);
                         cm = confusionmat(vallabel,str2num(cell2mat(predlabel)))%Compute the confusion matrix...
                         %Compute the accuracy mat...
                         accu_mat(kidx,1)=mean(vallabel==str2num(cell2mat(predlabel)));
                         disp(['Accuracy: ' num2str(accu_mat(kidx,1))]);
                         varimportance = ((kidx-1)*varimportance + mdl.OOBPermutedVarDeltaError)/kidx;
                         figure(1);
                         bar(varimportance);
                         title('Plot the variable importance...');
                         drawnow;
                        
                        
                 end
        end
        accu = mean(accu_mat);
        dev = std(accu_mat);
        disp(['Average accuracy: ' num2str(accu) ' +/- ' num2str(dev)]);
        
    end
end