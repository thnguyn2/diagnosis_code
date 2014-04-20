function [inmodel,history]=feature_selection_for_diagnosis(data,label)
%This function perform feature selections based on minimizing the
%classification error...

    opts = statset('display','iter','Display','iter');
    %Do sequential feature selection...
    c = cvpartition(size(data,1),'k',10);%Create cross-validation set
    misclasfunc = @(Xtrain,Ytrain,Xtest,Ytest)compute_misclassification_error(Xtrain,Ytrain,Xtest,Ytest);
    [inmodel,history]=sequentialfs(misclasfunc,data,label,'cv',c,'options',opts,'direction','backward');%The cross-validation is done here....
    %compute_misclassification_error(data,label,'knn');
end

function err=compute_misclassification_error(Xtrain,Ytrain,Xtest,Ytest)
%Compute the miss classification error based on k-fold validation
%Each row is an observation...
%First, create the k-fold validation set
    type = 'knn'; %Define the type of classifier...
    %Next, partion the training set
    switch type
        case {'knn'}
            %Build the k-NN classification
            mdl = ClassificationKNN.fit(Xtrain,Ytrain,'NumNeighbors',1);
            %Calculate the prediction error
            Ypred =predict(mdl,Xtest); %Calculate
            err = mean(Ypred==Ytest);           
        case {'tree'}
            %Build classifcation tree
            mdl = ClassificationTree.fit(Xtrain,Ytrain);
            %Calculate the prediction error
            Ypred =predict(mdl,Xtest); %Calculate
            err = mean(Ypred==Ytest);  
            
    end

end