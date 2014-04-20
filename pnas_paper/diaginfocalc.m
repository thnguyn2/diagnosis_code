clc;
clear all;
close all;
%This program compute the diagnosis information for different images in the
%dataset
disp('Calculating diagnosis information...');
labelpath = 'H:\QPI data\svm\texton_ws40_1024\';
fimpath = 'H:\QPI data\svm\texton_ws40_1024\';
morppath = 'H:\QPI data\morpinfo\glandsize\';
addpath(labelpath);
addpath(fimpath);
addpath('C:\Users\thnguyn2\Dropbox\current working code\Bilateral filtering'); %Add the path for bilateral filtering
datapath = 'H:\QPI data\'; %This is the path to all the SLIM images
textonpath = strcat(datapath,'texdir\');

[filenames,glandnames]=findFileName(datapath);

% %% ==================Classification based on the morphology of the tissue====================
load(strcat(morppath,'glandmorpbf.mat'));
featnorm = 1;
n33 =26;
n44 = 21;
nsamplestotal=67;
n33start = 1;
n44start = nsamplestotal-(n44-1);
n33train = round(n33*0.75);
n44train = round(n44*0.75);
n33test = n33 - n33train;
n44test = n44 - n44train;
%Normalize all the features to have zero-mean and std of 1
meanvect = mean(feat,1);
stdvect = std(feat,1);
nsamples = size(feat,1);
ndim = size(feat,2);

%Use the morphological feature
nsamples = size(feat,1);
ndim = size(feat,2);
if (featnorm==1) %If we do feature normalization
    meanvect = mean(feat,1);
    stdvect = std(feat,1);
    normfeat = (feat-repmat(meanvect,[nsamples 1]))./repmat(stdvect+1e-6,[nsamples 1]);
else
    normfeat = feat;
end
class33idx = find(dataclass==1);
class44idx = find(dataclass==4);

dataset = normfeat([class33idx(:);class44idx(:)],:);
datalabel = dataclass([class33idx(:);class44idx(:)],:);

%Randomly permute the data so that data from different classes are mixed...
p = randperm(length(datalabel));
dataset = dataset(p,:);
datalabel = datalabel(p,:);
disp('Classification using all features with knn')
[accu,dev,mdl,varimp] = estimate_kfold_prediction_error(dataset(:,:),datalabel,5,'knn');


% disp('Classification using all features with an SVM classifier')
% [accu,dev] = estimate_kfold_prediction_error(dataset,datalabel,5,'svm');


%Compute the k-fold prediction error when all features are used...
disp('Classification using all features with tree bagging')
[accu,dev,mdl,varimp] = estimate_kfold_prediction_error(dataset,datalabel,5,'treebagger');
[val,idx]=sort(varimp,'descend');
figure(5);
bar(varimp);
title('Variable importance...');
nmostimpt = 5;
%Plot the feature of two most important features..
feat1idx = idx(1);
feat2idx = idx(2);
figure(1);
for n33idx =1:n33
    plot(normfeat(class33idx(n33idx),feat1idx),normfeat(class33idx(n33idx),feat2idx),'+r');
    hold on;    
end
for n44idx =1:n44
    plot(normfeat(class44idx(n44idx),feat1idx),normfeat(class44idx(n44idx),feat2idx),'ob');
    hold on;    
end
title('Two most important features scatter plot');
xlabel(strcat('Feat ',num2str(feat1idx)));
ylabel(strcat('Feat ',num2str(feat2idx)));


%Perform classification with 10 most importance features...
disp('Classification with knn for 10 most importance features')
[accu,dev,mdl,varimp] = estimate_kfold_prediction_error(dataset(:,idx(1:nmostimpt)),datalabel,10,'knn');

%Next, perfrom feature selection and keep only a few features.
[infeat]=feature_selection_for_diagnosis(dataset,datalabel);
goodfeatidx = find(infeat==1);
%Estimate the error when some feature is removed to reduce predictor's
%variance
[accufs,devfs] = estimate_kfold_prediction_error(dataset(:,goodfeatidx),...
    datalabel,10,'knn');

