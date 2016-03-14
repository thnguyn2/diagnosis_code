function [conf, result] = multisvm(TrainingSet,GroupTrain,TestSet)
%This source code perform muti-class SVM classification based on 1 vs all
%criteria.
%Last edits:
%   Tan H. Nguyen
%   University of Illinois at Urbana-Champaign
%Intputs:
%   TrainingSet: a n x k matrix where each row is an observations. n observation totally.
%   k is the dimensions of the feature vector
%   Group train: a column vector where each row contains a label for each
%   observation. Totally, there are L group
%   Testset: a m x k matrix where each row is an observation for testing
%Outputs:
%   conf: a m x 1 matrix where each row is for one testing sample, each
%   column corresponds to a different classes
%   results: a m x 1 column vector that is result[m] = argmax_k(pmap(m,k))
%Notes:
%Models a given training set with a corresponding group vector and 
%classifies a given test set using an SVM classifier according to a 
%one vs. all relation. 
%
%Modified from the original code by Cody Neuburger cneuburg@fau.edu
%Florida Atlantic University, Florida USA
%This code was adapted and cleaned from Anand Mishra's multisvm function
%found at http://www.mathworks.com/matlabcentral/fileexchange/33170-multi-class-support-vector-machine/

u=unique(GroupTrain);
numClasses=length(u);
%build models
for k=1:numClasses
    %Vectorized statement that binarizes Group
    %where 1 is the current class and 0 is all other classes
    G1vAll=(GroupTrain==u(k));
    models(k) = svmtrain(TrainingSet,G1vAll,'kernel_function','rbf');
end
m = size(TestSet,1);
l = length(u(:));%Number of classes.
pmap = zeros(m,l);
%Go through all the classifier and compute the score of all the classes
for ll = 1:l 
        %Load the value of a classifier to compute
        shiftval = models(ll).ScaleData.shift;
        scaleval = models(ll).ScaleData.scaleFactor;
        supvect = models(ll).SupportVectors;%Get the support vectors
        alpha = models(ll).Alpha; %Note that alpha is positive for the 1st group and -1 for the second group
        bias = models(ll).Bias;
        kerfunc = models(ll).KernelFunction;
        kerfuncargs = models(ll).KernelFunctionArgs;
        curevalset =  TestSet;                       
        curevalset = bsxfun(@plus,curevalset,shiftval);%Data normalization
        curevalset = bsxfun(@times,curevalset,scaleval);
        pmap(:,ll)=-(kerfunc(supvect,curevalset,kerfuncargs{:})'*alpha(:) + bias(:));
                              
end

[conf,idx] = max(pmap,[],2);
result = u(idx);
end