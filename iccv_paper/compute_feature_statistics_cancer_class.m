function [meanclass2,meanclass3,meanclass4,meanclass5,stdclass2,stdclass3,stdclass4,stdclass5,cancergrade,class2feat,class3feat,class4feat,class5feat,c2name,c3name,c4name,c5name...
    ]=compute_feature_statistics_cancer_class(featvect,dataclass,donorm,samplename)
    %Compute the mean and the standard deviation of all features for 5
    %different groups.
    %feat: a 2D array where each row is a data sample, each column is a
    %feature vector.
    meanfeat = mean(featvect,1);
    if (nargin==2)
        donorm=false;
    end
    if (donorm==true)
            %Normalization every features so that they have zero mean and
            %std of 1
            featvect = featvect-repmat(meanfeat,[size(featvect,1) 1]);
            featvect = featvect./repmat(std(featvect,1),[size(featvect,1) 1]);
            %featvect = featvect./repmat(meanfeat,[size(featvect,1) 1]); %Normalize every feature so that every feature has mean value of 1
    end
    nfeat = size(featvect,2);
    meanclass2 = zeros(1,nfeat);
    meanclass3 = zeros(1,nfeat);
    meanclass4 = zeros(1,nfeat);
    meanclass5 = zeros(1,nfeat);
    stdclass2 = zeros(1,nfeat);
    stdclass3 = zeros(1,nfeat);
    stdclass4 = zeros(1,nfeat);
    stdclass5 = zeros(1,nfeat);
    class2feat = zeros(0,nfeat);
    class3feat = zeros(0,nfeat);
    class4feat = zeros(0,nfeat);
    class5feat = zeros(0,nfeat);
    c2name = zeros(1,0);
    c3name = zeros(1,0);
    c4name = zeros(1,0);
    c5name = zeros(1,0);
    cancergrade = zeros(size(dataclass));
    idxclass2=[  find(strcmp(dataclass,'22')==1);....
                 find(strcmp(dataclass,'23')==1)...
              ];  
    if (~isempty(idxclass2))
        meanclass2=mean(featvect(idxclass2,:),1);
        stdclass2 =std(featvect(idxclass2,:),0,1);
        cancergrade(idxclass2)=2;
        class2feat = featvect(idxclass2,:);
        c2name = samplename(idxclass2);
        
    end
    
    idxclass3 = [find(strcmp(dataclass,'33')==1);
                 find(strcmp(dataclass,'32')==1);
                 find(strcmp(dataclass,'34')==1)
              ];  
    if (~isempty(idxclass3))
        meanclass3=mean(featvect(idxclass3,:),1);
        stdclass3 =std(featvect(idxclass3,:),0,1);
        cancergrade(idxclass3)=3;
        class3feat = featvect(idxclass3,:);
        c3name = samplename(idxclass3);
    end
    
    idxclass4 = [find(strcmp(dataclass,'43')==1);
                 find(strcmp(dataclass,'44')==1);...
                 find(strcmp(dataclass,'45')==1)...
              ]; 
    if (~isempty(idxclass4))
        meanclass4=mean(featvect(idxclass4,:),1);
        stdclass4 = std(featvect(idxclass4,:),0,1);
        cancergrade(idxclass4)=4;
        class4feat = featvect(idxclass4,:);
        c4name = samplename(idxclass4);
    end
    
   
    idxclass5=[  find(strcmp(dataclass,'54')==1);....
                 find(strcmp(dataclass,'55')==1);...
                 find(strcmp(dataclass,'53')==1);....
                 find(strcmp(dataclass,'54')==1);....
                 find(strcmp(dataclass,'55')==1);];
    if (~isempty(idxclass5))
        meanclass5=mean(featvect(idxclass5,:),1);
        stdclass5 =std(featvect(idxclass5,:),0,1);
        cancergrade(idxclass5)=5;
        class5feat = featvect(idxclass5,:);
        c5name = samplename(idxclass5);
    end
   
end